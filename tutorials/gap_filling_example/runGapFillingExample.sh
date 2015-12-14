#assume you are inside the tutorials/gap_filling_example folder
#GapSeq.fa contains the gap sequence (Ns) and its flanking sequences (specifically, chr17:489695-491811 on hg38)
#PredictedGap.fa contains the sequence that is predicted to fill the above gap region (specifically, 000850F-001-01:2615-4744 on hx1 assembly)
../../bin/gfa.pl GapSeq.fa PredictedGap.fa > GapSeq.output
echo 'Please check GapSeq.output for results'
