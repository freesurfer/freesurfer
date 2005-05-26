LoadVolume "/home/kteich/test_data/anatomical/25cubed.mgz" 1 0
SetLayerLabel 0 "Layer 0"

for  { set n 1 } { $n < 10 } { incr n } {
    Make2DMRILayer "Layer $n"
    Set2DMRILayerVolumeCollection $n 0
    SetLayerInAllViewsInFrame 0 $n
}
