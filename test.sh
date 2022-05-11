parallel --dry-run -k -j 11 "./maxmeandp-vlns -i {1} -t 600 -s {2}" ::: inst/*.txt > t3-600s.log ::: `seq 1 5`
parallel --dry-run -k -j 11 "./mammdp_fv2 {1} 600 {2}" ::: inst/*.txt > t4-600s.log ::: `seq 1 5`

parallel --dry-run -k -j 11 "./maxmeandp-vlns -i {1} -t 1800 -s {2}" ::: inst/*.txt > t3-1800s.log ::: `seq 1 5`
parallel --dry-run -k -j 11 "./mammdp_fv2 {1} 1800 {2}" ::: inst/*.txt > t4-1800s.log ::: `seq 1 5`

