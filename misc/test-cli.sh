#!/bin/bash
set -x
sk_a="$(sidh csidh-genkey)"
pk_a="$(echo "$sk_a"|sidh csidh-pubkey -)"
sk_b="$(sidh csidh-genkey)"
pk_b="$(echo "$sk_b"|sidh csidh-pubkey -)"
ss_a="$(echo "$sk_a"|sidh csidh-dh - "$pk_b")"
ss_b="$(echo "$sk_b"|sidh csidh-dh - "$pk_a")"
if [ $? -eq 0 ] && [ "$ss_a" == "$ss_b" ]; then
    echo "OK"
else
    echo "fail: $sk_a $sk_b $pk_a $pk_b $ss_a $ss_b"
    exit 1
fi
