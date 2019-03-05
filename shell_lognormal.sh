rm OUT*
rm out*

cp sq_digestive.txt sq.txt

time ./lognormalh3>out_digestive<<@
0
-99
100000	N
20	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
0	EPIS
0.5	h_a
1.0	VE
2	REP
@
cat out_digestive >> OUT_DIGESTIVE

time ./lognormalh3>out_digestive<<@
0
-99
100000	N
20	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
1	EPIS
0.5	h_a
1.0	VE
2	REP
@
cat out_digestive >> OUT_DIGESTIVE

time ./lognormalh3>out_digestive<<@
0
-99
100000	N
20	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
0	EPIS
0.2	h_a
1.0	VE
2	REP
@
cat out_digestive >> OUT_DIGESTIVE

time ./lognormalh3>out_digestive<<@
0
-99
100000	N
20	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
1	EPIS
0.2	h_a
1.0	VE
@
cat out_digestive >> OUT_DIGESTIVE

