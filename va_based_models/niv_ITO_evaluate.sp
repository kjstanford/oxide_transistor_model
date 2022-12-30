* VSCNFET Demo: Id-Vds and Id-Vgs of N-type CNFET

.options POST

.param TEMP=273.13+27
.param L=1.5
.param W=20
.param R=500/W
.hdl 'ITO_TFT_compact_v5.va'

.param supply_max=2
.param supply_min=-0.5

xnfet1 d g s ITO_TFT_compact_v5 mu_band=39.5e-4 Rs=R Rd=R T=TEMP L=L W=W
* xnfet2 d g s ITO_TFT_compact_v5 mu_band=44.5e-4 Rs=R Rd=R T=TEMP L=L W=W

Vd d gnd 0.1
Vg g gnd 0
Vs s gnd 0

***********************************************************************
* Measurements
***********************************************************************
* test xnFET Ids vs. Vgs
.DC       Vg  START='supply_min'     STOP='supply_max'   STEP='0.02'

* .print dc idd1=par('i(xnfet1.d)') idd2=par('i(xnfet2.d)')
.print dc idd=par('i(xnfet1.d)')

.end