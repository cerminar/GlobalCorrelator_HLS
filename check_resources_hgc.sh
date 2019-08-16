#!/bin/sh

mkdir -p reports_hgc_alt

#clock=(5.000 4.167 3.571 3.125 2.500)
#clock_freq=(200.000 240.000 280.000 320.000 400.000)

#clock=(5.000 4.167 3.571 3.125)
#clock_freq=(200.000 240.000 280.000 320.000)

#clock=(5.000 3.571 3.125)
#clock_freq=(200.000 280.000 320.000)

clock=(4.167)
clock_freq=(240.000)

for ntrk in 38
do
  for ncalo in 38
  do
    for ii in 2
    do
      for iclk in 0
      do
        #for part in 'xcku115-flvb2104-2-i'
        for part in 'xcvu9p-flgb2104-2-i'
        do
          clockname=`echo "${clock_freq[iclk]}" | sed -r 's/\..*//g'`
          partname=`echo "$part" | sed -r 's/-.*//g'`
          echo $clockname
          echo $partname
          ./run_barebones_hgc.sh ${ntrk} ${ncalo} ${ii} ${clock[iclk]} ${clock_freq[iclk]} ${part}
          cp l1pfpuppi-hgc-resource-test/solution/syn/report/mp7wrapped_pfalgo2_hgc_csynth.rpt reports_hgc_alt/mp7wrapped_pfalgo2_hgc_csynth_ntrk${ntrk}_ncalo${ncalo}_ii${ii}_clock${clockname}_${partname}.rpt
          cp pftest_tieoff/pftest_tieoff.runs/impl_1/design_1_wrapper_power_routed.rpt reports_hgc_alt/design_1_wrapper_power_routed_ntrk${ntrk}_ncalo${ncalo}_ii${ii}_clock${clockname}_${partname}.rpt
          #cp puppi/proj_pfpuppi/solution1/syn/report/simple_puppi_hw_csynth.rpt reports_puppi/simple_puppi_hw_csynth_ntrk${ntrk}_ncalo${ncalo}.rpt
        done
      done
    done
  done
done
