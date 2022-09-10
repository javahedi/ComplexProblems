file="occup_dist_N10_2"


datafile1="Nph_2_dt_01_gs_occup_N10_mBt8.dat"
datafile2="Nph_2_dt_01_gs_occup_N10_mBt9.dat"
datafile3="Nph_2_dt_01_gs_occup_N10_mBt10.dat"
datafile4="Nph_2_dt_01_gs_occup_N10_mBt5.dat"
datafile5="Nph_2_dt_01_gs_occup_N10_mBt16.dat"
datafile6="Nph_2_dt_01_gs_occup_N10_mBt6.dat"
datafile8="Nph_2_dt_01_gs_occup_N20_mBt8.dat"
datafile10="Nph_2_dt_01_gs_occup_N20_mBt20.dat"
gnuplot << EOF

set term post eps enh color font ",16"
set output "$file.eps"
set title "Ground state particle occupation at t=20"
set xlabel "i"
set ylabel "n_i"
set key  right bottom
#set yr[.27:.3]	
       
   
     pl     "$datafile2"  u  (\$1+2):(\$2==200 ? \$3 :1/0)  w lp  t   "mBt=9",\
            "$datafile5"  u  (\$1+2):(\$2==200 ? \$3 :1/0)  w lp  t   "mBt=16",\
            "$datafile1"  u  (\$1+2):(\$2==200 ? \$3 :1/0)  w lp  t   "mBt=8",\
            "../../Dicke_exact_diagonalization/tight_binding/figures/occupation_dist_obc_M10.dat" u 1:202 w lp t   "exact"
           # "../../Dicke_exact_diagonalization/tight_binding/figures/occupation_dist_obc2.dat" u 1:102 w lp t   "exact2"
      
EOF

epstopdf $file.eps
rm -f $file.eps
