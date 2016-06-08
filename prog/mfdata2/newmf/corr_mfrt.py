#!/usr/bin/env python

import sys, os, re

fnin = "mfrt.log"

s = open(fnin).readlines()
sout = ""
for ln in s:
  m = re.match(r"(.*)\| torq (.+) (.+) (.+)\| symmtorq (.+) (.+) (.+)", ln.strip())
  if not m:
    print ln
    raise Exception
  symmtorq = float( m.group(2) ) * 0.5
  symmtorq_std = float( m.group(3) ) * 0.5
  torq = float( m.group(5) ) * 2
  torq_std = float( m.group(6) ) * 2
  sout += "%s| torq %s %s %s| symmtorq %s %s %s\n" % (
      m.group(1), torq, torq_std, m.group(7),
      symmtorq, symmtorq_std, m.group(4) )

print sout
