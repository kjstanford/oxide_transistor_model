* VSCNFET Demo: Id-Vds and Id-Vgs of N-type CNFET

.options POST

.param TEMP=25
.hdl 'ITO_TFT_compact_v5.va'

.param supply_max=2
.param supply_min=-0.5

xnfet d g s ITO_TFT_compact_v5

Vd d gnd 0.1
Vg g gnd 0
Vs s gnd 0

***********************************************************************
* Measurements
***********************************************************************
* test xnFET Ids vs. Vgs
.DC       Vg  START='supply_min'     STOP='supply_max'   STEP='0.01'

.print dc idd=par('i(xnfet.d)')

.end