SetFrameViewConfiguration 1 c1
set colID [MakeDataCollection Volume]
SetVolumeCollectionFileName $colID /home/kteich/freesurfer/subjects/bert/mri/T1
set layerID [MakeLayer 2DMRI]
SetLayerLabel $layerID "bert"
Set2DMRILayerVolumeCollection $layerID $colID
set viewID [GetViewIDFromFrameColRow 1 0 0]
set level [GetFirstUnusedDrawLevelInView $viewID]
SetLayerInViewAtLevel $viewID $layerID $level
