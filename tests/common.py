# -*- coding: utf-8 -*-
# FloodDepth - Licensed under the GNU General Public License v3 or later.
"""Tiny check harness shared by the pure suites (QEHT convention:
each check prints [PASS] / [FAIL]; the suite exits 1 on any failure)."""

import sys

FAILURES = []


def check(name, condition, detail=""):
    print(f"  [{'PASS' if condition else 'FAIL'}] {name}" + (f"   {detail}" if detail else ""))
    if not condition:
        FAILURES.append(name)


def finish(suite):
    print()
    if FAILURES:
        print(f"{suite}: {len(FAILURES)} FAILED: {FAILURES}")
        sys.exit(1)
    print(f"{suite}: all checks passed")
