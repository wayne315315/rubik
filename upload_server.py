"""Phone-friendly upload page that runs the photo pipeline on wayne-kv.

Open http://<wayne-kv>:8080 on the phone (LAN or NordVPN meshnet), take or pick
the two corner-view photos, tap Upload. The photos are saved under
uploads/<timestamp>/, cube_vision runs on them with the default Ollama model,
and the page shows the answer and scramble sequences plus a link to answer.mp4.
Standard library only; run with the project's virtualenv:

    python upload_server.py            # 0.0.0.0:8080
    python upload_server.py --port 9000 --optimal
"""
import argparse
import html
import json
import os
import subprocess
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from email.parser import BytesParser
from email.policy import default as email_policy

ROOT = os.path.dirname(os.path.abspath(__file__))
UPLOADS = os.path.join(ROOT, "uploads")
JOBS = {}                                    # job id -> dict(status, log, result)
EXTRA_ARGS = []

PAGE = """<!doctype html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Cube upload</title>
<style>body{font-family:-apple-system,sans-serif;margin:1.2em;max-width:40em}input,button{font-size:1.1em;margin:.4em 0}
pre{white-space:pre-wrap;background:#f4f4f4;padding:.6em;border-radius:.4em}</style>
<h2>Cube photos</h2>
<form method="post" enctype="multipart/form-data">
<p>Photo 1 (three faces): <input type="file" name="p1" accept="image/*" required></p>
<p>Photo 2 (the other three faces): <input type="file" name="p2" accept="image/*" required></p>
<button type="submit">Upload and solve</button></form>
<p>Recent jobs:</p><ul>%s</ul>"""

JOB_PAGE = """<!doctype html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
%s<title>Job %s</title>
<style>body{font-family:-apple-system,sans-serif;margin:1.2em;max-width:40em}pre{white-space:pre-wrap;background:#f4f4f4;padding:.6em;border-radius:.4em}</style>
<h2>Job %s: %s</h2>%s<pre>%s</pre><p><a href="/">Back</a></p>"""


def run_job(job_id, folder):
    job = JOBS[job_id]
    cmd = [sys.executable, "-u", os.path.join(ROOT, "cube_vision.py"),
           os.path.join(folder, "photo1.jpg"), os.path.join(folder, "photo2.jpg"),
           "--save-json", os.path.join(folder, "reading.json"), "--video"] + EXTRA_ARGS
    env = dict(os.environ, MPLBACKEND="Agg")
    proc = subprocess.Popen(cmd, cwd=folder, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
    with open(os.path.join(folder, "log.txt"), "w") as fh:          # keep the log on disk too
        for line in proc.stdout:
            job["log"] += line
            fh.write(line)
            fh.flush()
    proc.wait()
    job["status"] = "done" if proc.returncode == 0 else "failed"
    for line in job["log"].splitlines():
        if line.startswith(("Answer", "Scramble")):
            job["result"].append(line)


class Handler(BaseHTTPRequestHandler):
    def _send(self, body, ctype="text/html; charset=utf-8", code=200):
        data = body if isinstance(body, bytes) else body.encode()
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if self.path == "/":
            items = "".join(f'<li><a href="/job/{j}">{j}</a> {JOBS[j]["status"]}</li>' for j in sorted(JOBS, reverse=True)[:20])
            return self._send(PAGE % (items or "<li>none yet</li>"))
        if self.path.startswith("/job/"):
            job_id = self.path.split("/")[2]
            job = JOBS.get(job_id)
            if not job:
                return self._send("no such job", code=404)
            refresh = '<meta http-equiv="refresh" content="5">' if job["status"] == "running" else ""
            video = f'<p><video controls playsinline src="/video/{job_id}" style="max-width:100%%"></video></p>' if job["status"] == "done" else ""
            result = "<br>".join(html.escape(r) for r in job["result"])
            return self._send(JOB_PAGE % (refresh, job_id, job_id, job["status"],
                                          f"<p><b>{result}</b></p>{video}" if result else video, html.escape(job["log"])))
        if self.path.startswith("/video/"):
            job_id = self.path.split("/")[2]
            path = os.path.join(UPLOADS, job_id, "answer.mp4")
            if job_id in JOBS and os.path.exists(path):
                with open(path, "rb") as fh:
                    return self._send(fh.read(), "video/mp4")
        self._send("not found", code=404)

    def do_POST(self):
        length = int(self.headers.get("Content-Length", 0))
        body = self.rfile.read(length)
        msg = BytesParser(policy=email_policy).parsebytes(
            b"Content-Type: " + self.headers["Content-Type"].encode() + b"\r\n\r\n" + body)
        files = {part.get_param("name", header="content-disposition"): part.get_payload(decode=True)
                 for part in msg.iter_parts()}
        if not files.get("p1") or not files.get("p2"):
            return self._send("need two photos", code=400)
        job_id = time.strftime("%Y%m%d-%H%M%S")
        folder = os.path.join(UPLOADS, job_id)
        os.makedirs(folder, exist_ok=True)
        for name, key in (("photo1.jpg", "p1"), ("photo2.jpg", "p2")):
            with open(os.path.join(folder, name), "wb") as fh:
                fh.write(files[key])
        JOBS[job_id] = {"status": "running", "log": "", "result": []}
        threading.Thread(target=run_job, args=(job_id, folder), daemon=True).start()
        self.send_response(303)
        self.send_header("Location", f"/job/{job_id}")
        self.end_headers()

    def log_message(self, fmt, *args):
        sys.stderr.write("%s - %s\n" % (self.address_string(), fmt % args))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=8080)
    parser.add_argument("--optimal", action="store_true", help="pass --optimal to cube_vision")
    parser.add_argument("--model", default=None, help="model name to pass to cube_vision")
    args = parser.parse_args()
    if args.optimal:
        EXTRA_ARGS.append("--optimal")
    if args.model:
        EXTRA_ARGS.extend(["--model", args.model])
    os.makedirs(UPLOADS, exist_ok=True)
    print(f"serving on http://{args.host}:{args.port}  (uploads -> {UPLOADS})")
    ThreadingHTTPServer((args.host, args.port), Handler).serve_forever()


if __name__ == "__main__":
    main()
