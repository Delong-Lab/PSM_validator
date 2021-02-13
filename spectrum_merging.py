'''
run the new data (and can try with old too) using msconvert-converted vs specmill-converted; how do the results differ?
if behaves like you expect, figure out what specmill does that helps
is it precursor refine? if so, see how differently precursor refine works for the two programs; do i need to write my own precursor refine? i think you would just pick the peak closest in mass and that is biggest
why doesn't vendor software wait and do the mass accuracy calibration after everything is done so it could easily know which parent peak the ms2 metadata precursor corresponds to and just change it with that?
It's because the quad tried to grab that spot, so you basically just want to pick whatever parent ion falls closest to where the quad targeted!!!!!!!!!!!!!!!!
is it merging? if so, write your own merging method in which you sum intensities and average the m/zs; could use PCC_calculator as basis for this
you would also basically keep a running list of all precursors with m/z and RT. 
with each new one, you look back at the list and once you are in the right RT range see if any of the previous scans have the same m/z
if they do, could then check for some level of similarity and if they pass that then you merge the scans
NOTE: the analyses seem to work just fine on the benchmarking data, so what is different about that data and the hip data that doesn't work? i think it's the abundance of the bio hip maybe
'''