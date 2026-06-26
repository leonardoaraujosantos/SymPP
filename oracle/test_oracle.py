#!/usr/bin/env python3
"""Regression tests for the SymPP oracle harness.

Run with:  python3 -m unittest discover -s oracle -p 'test_*.py'
       or:  just test-oracle-py
"""
import json
import os
import subprocess
import sys
import unittest

ORACLE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "oracle.py")


class ShutdownTest(unittest.TestCase):
    def test_shutdown_with_closed_reader_does_not_traceback(self):
        """The C++ harness sends "shutdown" then closes the read end of the
        pipe without consuming the reply. The oracle must exit cleanly instead
        of dying with a BrokenPipeError traceback on stderr.
        """
        proc = subprocess.Popen(
            [sys.executable, ORACLE],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        # Mirror tests/oracle/oracle.cpp::shutdown_subprocess(): write the
        # shutdown request, then drop the read end straight away.
        proc.stdin.write((json.dumps({"op": "shutdown"}) + "\n").encode())
        proc.stdin.flush()
        proc.stdout.close()  # parent stops reading replies
        proc.stdin.close()
        self.assertEqual(proc.wait(timeout=15), 0, f"oracle exited {proc.returncode}")
        stderr = proc.stderr.read()
        proc.stderr.close()

        self.assertNotIn(b"Traceback", stderr, stderr.decode(errors="replace"))
        self.assertNotIn(b"BrokenPipeError", stderr, stderr.decode(errors="replace"))

    def test_ping_then_shutdown(self):
        """A normal request still gets a reply before shutdown ends the loop."""
        proc = subprocess.Popen(
            [sys.executable, ORACLE],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, stderr = proc.communicate(
            (json.dumps({"id": 1, "op": "ping"}) + "\n"
             + json.dumps({"id": 2, "op": "shutdown"}) + "\n").encode(),
            timeout=15,
        )
        lines = [l for l in out.decode().splitlines() if l.strip()]
        self.assertEqual(len(lines), 1, f"expected only the ping reply, got {lines!r}")
        reply = json.loads(lines[0])
        self.assertTrue(reply["ok"])
        self.assertEqual(reply["result"], "pong")
        self.assertEqual(proc.returncode, 0)
        self.assertNotIn(b"Traceback", stderr)


if __name__ == "__main__":
    unittest.main()
