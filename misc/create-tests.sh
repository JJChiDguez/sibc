#!/bin/bash

# this generates known value tests for various parameters

create_known_answer_script() {
    echo "echo 'starting test of $@'"

    sk_a="$(sidh "$@" genkey)"
    echo "sk_a='$sk_a'"
    pk_a="$(echo "$sk_a"|sidh "$@" pubkey -)"
    echo "pk_a='$pk_a'"
    echo '[ "$pk_a" != "$(echo "$sk_a"|sidh '"$@"' pubkey -)" ] && fail=1 && echo "fail: pk($sk_a) -> $pk_a" || echo "OK: pk_a '"$@"'"'

    sk_b="$(sidh "$@" genkey)"
    echo "sk_b='$sk_b'"
    pk_b="$(echo "$sk_b"|sidh "$@" pubkey -)"
    echo "pk_b='$pk_b'"
    echo '[ "$pk_b" != "$(echo "$sk_b"|sidh '"$@"' pubkey -)" ] && fail=1 && echo "fail: pk($sk_b) -> $pk_b" || echo "OK: pk_b '"$@"'"'

    ss_a="$(echo "$sk_a"|sidh "$@" dh - "$pk_b")"
    echo "ss_a='$ss_a'"
    echo '[ "$ss_a" != "$(echo "$sk_a"|sidh '"$@"' dh - "$pk_b")" ] && fail=1 && echo "fail: dh($sk_a,$pk_b) -> $ss_a" || echo "OK: ss_a '"$@"'"'

    ss_b="$(echo "$sk_b"|sidh "$@" dh - "$pk_a")"
    echo "ss_b='$ss_b'"
    echo '[ "$ss_b" != "$(echo "$sk_b"|sidh '"$@"' dh - "$pk_a")" ] && fail=1 && echo "fail: dh($sk_b,$pk_a) -> $ss_b" || echo "OK: ss_b '"$@"'"'

    [ "$ss_a" != "$ss_b" ] && echo "fail: $@" && exit 1

    echo '[ "$ss_a" != "$ss_b" ] && fail=1 && echo "fail: $ss_a != $ss_b" || echo "OK: same '"$@"'"'
    echo
}

echo '# script generated by create-tests.sh'
echo

for prime in p512 p1024 p1792; do
    echo
    echo "# $(date) # $prime"
    for formula in hvelu tvelu svelu; do
        for style in df wd1 wd2; do
            for tuned in "--tuned" ""; do
                for multievaluation in "--multievaluation" ""; do
                    if [ "$formula" == "tvelu" ] && ([ "$multievaluation" != "" ] || [ "$tuned" != "" ])
                    then
                        continue
                    fi
                    echo
                    create_known_answer_script -p $prime -f $formula -s $style $tuned $multievaluation
                done
            done
        done
    done
done
echo
echo "# create-tests.sh finished at $(date)"

