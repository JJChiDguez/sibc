#!/bin/bash
echo ""
echo "CSIDH"
set -x
sk_a="$(sibc csidh-genkey)"
pk_a="$(echo "$sk_a"|sibc csidh-pubkey -)"
sk_b="$(sibc csidh-genkey)"
pk_b="$(echo "$sk_b"|sibc csidh-pubkey -)"
ss_a="$(echo "$sk_a"|sibc csidh-dh - "$pk_b")"
ss_b="$(echo "$sk_b"|sibc csidh-dh - "$pk_a")"
if [ $? -eq 0 ] && [ "$ss_a" == "$ss_b" ]; then
    echo "OK"
else
    echo "fail: $sk_a $sk_b $pk_a $pk_b $ss_a $ss_b"
    exit 1
fi
echo ""
echo "SIKE"
sk="$(sibc -a sidh -p p434 keygen)"
pk3=`echo "${sk}" | tail -n1`
ck="$(echo "$pk3"|sibc -a sidh -p p434 encaps -)"
c0=`echo "${ck}" | head -1`
c1=`echo "${ck}" | tail -2 | head -1`
K=`echo "${ck}" | tail -n1`
K_="$(echo "$sk"|sibc -a sidh -p p434 decaps - "$c0 $c1")"
if [ $? -eq 0 ] && [ "$K" == "$K_" ]; then
    echo "OK"
else
    echo "fail: $sk $ck $K $K_"
    exit 1
fi
echo ""
echo "BSIKE"
sk="$(sibc -a bsidh -p p253 keygen)"
pk3=`echo "${sk}" | tail -n1`
ck="$(echo "$pk3"|sibc -a bsidh -p p253 encaps -)"
c0=`echo "${ck}" | head -1`
c1=`echo "${ck}" | tail -2 | head -1`
K=`echo "${ck}" | tail -n1`
K_="$(echo "$sk"|sibc -a bsidh -p p253 decaps - "$c0 $c1")"
if [ $? -eq 0 ] && [ "$K" == "$K_" ]; then
    echo "OK"
else
    echo "fail: $sk $ck $K $K_"
    exit 1
fi