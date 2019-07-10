#setenv XLSMPOPTS "parthds=6"
#setenv XLSMPOPTS "parthds=${1}:schedule=guided:spins=1:delays=1:yields=1"
time ../../main/odea <<EOF 
gen_hjs.sh
paramhjs.in
plhjs.in
tphjs.in
100
EOF
#
exit
