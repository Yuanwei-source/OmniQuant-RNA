#!/usr/bin/env python3
"""Shared helpers for parsing and normalizing annotation attributes."""

import re


ATTRIBUTE_SPLIT_RE = re.compile(r"\s*;\s*")


def parse_attributes(raw):
    """Parse GTF/GFF-style attribute field into a dict.

    Supports both key=value and key "value" styles.
    """
    attrs = {}
    for field in ATTRIBUTE_SPLIT_RE.split(str(raw).strip().strip(";")):
        token = field.strip()
        if not token:
            continue

        if "=" in token:
            key, value = token.split("=", 1)
        elif " " in token:
            key, value = token.split(" ", 1)
        else:
            continue

        attrs[key.strip()] = value.strip().strip('"')

    return attrs
