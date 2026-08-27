#!/bin/bash
declare -A LEN=( [1]=249 [2]=243 [3]=199 [4]=191 [5]=182 [6]=171 [7]=160 [8]=146 [9]=139 [10]=134
                 [11]=136 [12]=134 [13]=115 [14]=108 [15]=102 [16]=91 [17]=84 [18]=81 [19]=59
                 [20]=65 [21]=47 [22]=51 [X]=157 )
CHUNK=40
for c in "${!LEN[@]}"; do
  n=$(( (${LEN[$c]} + CHUNK - 1) / CHUNK ))
  for i in $(seq 0 $((n-1))); do
    s=$(( i*CHUNK*1000000 + 1 )); e=$(( (i+1)*CHUNK*1000000 ))
    f="c${c}_${i}.fa"
    [ -s "$f" ] && grep -q '^>' "$f" && continue
    Q="<?xml version=\"1.0\" encoding=\"UTF-8\"?><!DOCTYPE Query><Query virtualSchemaName=\"default\" formatter=\"FASTA\" header=\"0\" uniqueRows=\"1\" count=\"\" datasetConfigVersion=\"0.6\"><Dataset name=\"hsapiens_gene_ensembl\" interface=\"default\"><Filter name=\"chromosome_name\" value=\"${c}\"/><Filter name=\"start\" value=\"${s}\"/><Filter name=\"end\" value=\"${e}\"/><Filter name=\"transcript_is_canonical\" excluded=\"0\"/><Filter name=\"biotype\" value=\"protein_coding\"/><Attribute name=\"3utr\"/><Attribute name=\"external_gene_name\"/><Attribute name=\"ensembl_transcript_id\"/></Dataset></Query>"
    for try in 1 2 3 4 5; do
      curl -sSL -m 300 -o "$f" -G --data-urlencode "query=$Q" "https://www.ensembl.org/biomart/martservice"
      if grep -q '^>' "$f"; then echo "OK $f $(grep -c '^>' $f)"; break; fi
      sleep $((try*8))
    done
    grep -q '^>' "$f" || echo "FAIL $f"
    sleep 3
  done
done
echo "TOTAL $(cat c*.fa 2>/dev/null | grep -c '^>')"
