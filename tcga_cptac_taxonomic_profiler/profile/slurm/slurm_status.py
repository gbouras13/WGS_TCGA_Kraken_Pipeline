#!/usr/bin/env python3
"""Cluster status script for snakemake's --cluster-status.

Prints exactly one of: success | failed | running

WHY THIS IS NOT A ONE-LINER
---------------------------
The previous version asked sacct alone and treated anything it did not
recognise as a failure:

    output = sacct -j <id> --format State --noheader | head -1
    if "COMPLETED" in output:   success
    elif <a running state>:     running
    else:                       failed        # <-- including empty output

sacct reads SLURM's accounting database, which lags the controller by a few
seconds after submission. During that window sacct returns NOTHING for a
perfectly healthy job, the empty string matches no known state, and the job is
declared failed before it has even started.

That is exactly what happened on 2026-09-04. A COMEBin driver submitted its
child, polled faster than the accounting database updated, declared the job
failed, restarted it, declared that failed too, and exited 49 seconds later.
Both children then started normally three minutes afterwards and died only
because the driver had already torn down .snakemake/tmp.*. COMEBin itself ran
fine when invoked by hand. Nothing was wrong with the tool, the environment, or
the rule.

The risk scales with the number of jobs submitted at once, so it is worst
exactly where it hurts most: the full 336-sample run.

TWO CHANGES
-----------
1. Ask squeue first. squeue reads the controller, which knows about a job the
   instant it is submitted, and only fall back to sacct (needed for jobs that
   have already left the queue).
2. Never infer failure from absence. If neither source knows the job, retry
   briefly, then report "running" so snakemake polls again. A job wrongly
   reported as running is recoverable; a job wrongly reported as failed takes
   the whole workflow down with it.
"""

import shlex
import subprocess
import sys
import time

# Transient states: the job exists and has not finished.
RUNNING = {
    "PENDING", "CONFIGURING", "RUNNING", "COMPLETING", "SUSPENDED",
    "PREEMPTED", "RESIZING", "REQUEUED", "REQUEUE_HOLD", "REQUEUE_FED",
    "RESV_DEL_HOLD", "SIGNALING", "STAGE_OUT", "STOPPED",
}
SUCCESS = {"COMPLETED"}

# How long to keep asking before assuming accounting lag rather than failure.
ATTEMPTS = 5
BACKOFF = (1, 2, 3, 4)


def _capture(cmd):
    """Run cmd, return stripped stdout, or "" on any failure."""
    try:
        p = subprocess.run(
            shlex.split(cmd), capture_output=True, text=True, timeout=60
        )
    except Exception:
        return ""
    return p.stdout.strip() if p.returncode == 0 else ""


def _first_state(text):
    """First non-empty state token. Handles 'CANCELLED by 12345'."""
    for line in text.splitlines():
        tok = line.strip().split("|")[0].strip()
        if tok:
            return tok.split()[0].upper()
    return ""


def query(jobid):
    # Controller first - authoritative the moment a job is submitted.
    state = _first_state(_capture("squeue -j %s -h -o %%T" % jobid))
    if state:
        return state
    # Then accounting, for jobs that have already left the queue.
    return _first_state(
        _capture("sacct -j %s --format State --noheader --parsable2" % jobid)
    )


def main():
    if len(sys.argv) < 2:
        print("running")
        return
    jobid = sys.argv[-1].split(";")[0].strip()

    state = ""
    for attempt in range(ATTEMPTS):
        state = query(jobid)
        if state:
            break
        if attempt < len(BACKOFF):
            time.sleep(BACKOFF[attempt])

    if not state:
        # Neither source knows this job. Almost always accounting lag on a
        # freshly submitted job. Do NOT report failure from absence.
        print("running")
    elif state in SUCCESS:
        print("success")
    elif state in RUNNING:
        print("running")
    else:
        print("failed")


if __name__ == "__main__":
    main()
