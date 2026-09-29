# Digital-twin table check (v1822)

## digital_twin_properties_bcm_coreg_v4
schema: digital_twin_properties

description: 
Functional properties computed from the digital twin model from Wang et al. 2025, merged with coregistration_manual_v4.

    cc_abs : Test set performance of the digital twin model unit, higher is better
    cc_max : Neuron variability score used to normalize digital twin model unit performance
    cc_norm : Normalized model unit performance, higher is better
    OSI : orientation selectivity index
    DSI : direction selectivity index
    gOSI : global orientation selectivity index
    gDSI : global direction selectivity index
    pref_ori : Preferred orientation in degrees (0 - 180), vertical bar moving right is 0 and orientation increases counter-clockwise
    pref_dir : Preferred direction in degrees (0 - 360), vertical bar moving right is 0 and orientation increases counter-clockwise
    readout_loc_x : X coordinate of the readout location, an approximation of receptive field center in stimulus space; (-1, -1) top-left, (1, 1) bottom-right
    readout_loc_y : Y coordinate of the readout location, an approximation of receptive field center in stimulus space; (-1, -1) top-left, (1, 1) bottom-right

This table replaces `functional_properties_v3_bcm`
This table is maintained by Stelios Papadopoulos @ Tolias Lab

For the use of this table please cite:
1) https://doi.org/10.1038/s41586-025-08829-y; Wang et al., 2025
2) https://doi.org/10.1038/s41586-025-08840-3; Ding, Fahey, and Papadopoulos et al., 2025
3) https://doi.org/10.1038/s41586-025-08790-w; MICrONS Consortium et al., 2025
 [Note: This table 'digital_twin_properties_bcm_coreg_v4' will update the 'target_id' foreign_key when updates are made to the 'nucleus_detection_v0' table] 

columns: ['id_ref', 'created_ref', 'valid_ref', 'volume', 'pt_supervoxel_id', 'pt_root_id', 'id', 'created', 'valid', 'target_id', 'session', 'scan_idx', 'unit_id', 'pref_ori', 'pref_dir', 'gOSI', 'gDSI', 'cc_abs', 'cc_max', 'cc_norm', 'OSI', 'DSI', 'readout_loc_x', 'readout_loc_y', 'pt_position', 'bb_start_position', 'bb_end_position']

|    |   id_ref | created_ref                      | valid_ref   |   volume |   pt_supervoxel_id |         pt_root_id |   id | created                          | valid   |   target_id |   session |   scan_idx |   unit_id |   pref_ori |   pref_dir |      gOSI |      gDSI |    cc_abs |   cc_max |   cc_norm |      OSI |       DSI |   readout_loc_x |   readout_loc_y | pt_position            | bb_start_position   | bb_end_position   |
|---:|---------:|:---------------------------------|:------------|---------:|-------------------:|-------------------:|-----:|:---------------------------------|:--------|------------:|----------:|-----------:|----------:|-----------:|-----------:|----------:|----------:|----------:|---------:|----------:|---------:|----------:|----------------:|----------------:|:-----------------------|:--------------------|:------------------|
|  0 |   516506 | 2020-09-28 22:44:43.650751+00:00 | True        |  276.767 | 105487283075464806 | 864691135348268503 |    1 | 2025-07-30 17:51:22.691239+00:00 | True    |      516506 |         4 |          7 |       648 |    73.8787 |   253.879  | 0.104551  | 0.130178  | 0.357543  | 0.591973 |  0.603985 | 0.198416 | 0.379415  |      -0.14646   |      -0.273172  | [296208  94688  15803] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  1 |   452329 | 2020-09-28 22:45:02.852190+00:00 | True        |  311.618 | 101758289725398433 | 864691135700505634 |    2 | 2025-07-30 17:51:22.692288+00:00 | True    |      452329 |         4 |          7 |       662 |    92.6248 |   272.625  | 0.44016   | 0.0836913 | 0.470725  | 0.73643  |  0.639199 | 0.747158 | 0.182383  |      -0.0461872 |      -0.105799  | [269248  98832  18382] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  2 |   451461 | 2020-09-28 22:41:51.543636+00:00 | True        |  263.504 | 102531727972419182 | 864691135776919981 |    3 | 2025-07-30 17:51:22.693214+00:00 | True    |      451461 |         4 |          7 |       665 |    26.1577 |   206.158  | 0.0862164 | 0.0607932 | 0.290434  | 0.511465 |  0.567847 | 0.17052  | 0.252502  |      -0.122566  |      -0.185846  | [274864  94064  22046] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  3 |   420222 | 2020-09-28 22:45:01.445495+00:00 | True        |  309.883 | 100491446037118564 | 864691135472842290 |    4 | 2025-07-30 17:51:22.694580+00:00 | True    |      420222 |         4 |          7 |       671 |   111.083  |   291.083  | 0.387291  | 0.0284634 | 0.0696444 | 0.574129 |  0.121304 | 0.632664 | 0.0990792 |      -0.0904534 |      -0.0652829 | [259936  97360  16974] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  4 |   420058 | 2020-09-28 22:44:36.438460+00:00 | True        |  267.026 |  98661652329244694 | 864691135349237975 |    5 | 2025-07-30 17:51:22.695479+00:00 | True    |      420058 |         4 |          7 |       682 |    15.6009 |    15.6009 | 0.152346  | 0.0391141 | 0.324326  | 0.677929 |  0.478407 | 0.293705 | 0.141291  |      -0.0961552 |      -0.0826342 | [246448  95856  15873] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |

## digital_twin_properties_bcm_coreg_auto_phase3_fwd_v2
schema: digital_twin_properties

description: 
Functional properties computed from the digital twin model from Wang et al. 2025, merged with coregistration_auto_phase3_fwd_v2.

    cc_abs : Test set performance of the digital twin model unit, higher is better
    cc_max : Neuron variability score used to normalize digital twin model unit performance
    cc_norm : Normalized model unit performance, higher is better
    OSI : orientation selectivity index
    DSI : direction selectivity index
    gOSI : global orientation selectivity index
    gDSI : global direction selectivity index
    pref_ori : Preferred orientation in degrees (0 - 180), vertical bar moving right is 0 and orientation increases counter-clockwise
    pref_dir : Preferred direction in degrees (0 - 360), vertical bar moving right is 0 and orientation increases counter-clockwise
    readout_loc_x : X coordinate of the readout location, an approximation of receptive field center in stimulus space; (-1, -1) top-left, (1, 1) bottom-right
    readout_loc_y : Y coordinate of the readout location, an approximation of receptive field center in stimulus space; (-1, -1) top-left, (1, 1) bottom-right

This table is maintained by Stelios Papadopoulos @ Tolias Lab

For the use of this table please cite:
1) https://doi.org/10.1038/s41586-025-08829-y; Wang et al., 2025
2) https://doi.org/10.1038/s41586-025-08840-3; Ding, Fahey, and Papadopoulos et al., 2025
3) https://doi.org/10.1038/s41586-025-08790-w; MICrONS Consortium et al., 2025
 [Note: This table 'digital_twin_properties_bcm_coreg_auto_phase3_fwd_v2' will update the 'target_id' foreign_key when updates are made to the 'nucleus_detection_v0' table] 

columns: ['id_ref', 'created_ref', 'valid_ref', 'volume', 'pt_supervoxel_id', 'pt_root_id', 'id', 'created', 'valid', 'target_id', 'session', 'scan_idx', 'unit_id', 'pref_ori', 'pref_dir', 'gOSI', 'gDSI', 'cc_abs', 'cc_max', 'cc_norm', 'OSI', 'DSI', 'readout_loc_x', 'readout_loc_y', 'pt_position', 'bb_start_position', 'bb_end_position']

|    |   id_ref | created_ref                      | valid_ref   |   volume |   pt_supervoxel_id |         pt_root_id |    id | created                          | valid   |   target_id |   session |   scan_idx |   unit_id |   pref_ori |   pref_dir |      gOSI |      gDSI |   cc_abs |   cc_max |   cc_norm |       OSI |       DSI |   readout_loc_x |   readout_loc_y | pt_position            | bb_start_position   | bb_end_position   |
|---:|---------:|:---------------------------------|:------------|---------:|-------------------:|-------------------:|------:|:---------------------------------|:--------|------------:|----------:|-----------:|----------:|-----------:|-----------:|----------:|----------:|---------:|---------:|----------:|----------:|----------:|----------------:|----------------:|:-----------------------|:--------------------|:------------------|
|  0 |   106242 | 2020-09-28 22:41:12.393451+00:00 | True        |  261.994 |  79525545733567074 | 864691136309745498 | 35852 | 2025-07-30 18:08:20.155743+00:00 | True    |      106242 |         9 |          4 |      2388 |   124.556  |    304.556 | 0.210832  | 0.0323181 | 0.219483 | 0.59398  |  0.369512 | 0.40676   | 0.0263991 |       -0.214859 |       -0.20811  | [107104 127040  15087] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  1 |   106431 | 2020-09-28 22:43:14.353734+00:00 | True        |  289.367 |  79947689680746397 | 864691135725228203 | 35980 | 2025-07-30 18:08:20.253976+00:00 | True    |      106431 |         9 |          4 |      2536 |    65.9691 |    245.969 | 0.0231131 | 0.0842573 | 0.4073   | 0.504181 |  0.807845 | 0.0286045 | 0.282779  |       -0.252097 |       -0.143649 | [110288 126400  16890] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  2 |   106436 | 2020-09-28 22:41:20.272045+00:00 | True        |  315.426 |  80018127211293829 | 864691135133547168 | 35711 | 2025-07-30 18:08:20.056047+00:00 | True    |      106436 |         9 |          4 |      2212 |    20.3977 |    200.398 | 0.109068  | 0.0973767 | 0.482162 | 0.719846 |  0.669813 | 0.213756  | 0.276174  |       -0.202574 |       -0.196255 | [110768 126704  17166] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  3 |   106518 | 2020-09-28 22:41:14.258730+00:00 | True        |  309.334 |  80299602255318683 | 864691135210056384 | 35967 | 2025-07-30 18:08:20.244699+00:00 | True    |      106518 |         9 |          4 |      2521 |   103.407  |    283.407 | 0.118369  | 0.0749642 | 0.558106 | 0.830425 |  0.672073 | 0.231545  | 0.257083  |       -0.244328 |       -0.163041 | [112960 126672  17863] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  4 |   109045 | 2020-09-28 22:42:31.399120+00:00 | True        |  305.216 |  80440889163937016 | 864691136903060914 | 36818 | 2025-07-30 18:08:20.830762+00:00 | True    |      109045 |         9 |          4 |      5010 |   175.938  |    175.938 | 0.0339695 | 0.0619844 | 0.391196 | 0.624106 |  0.62681  | 0.104322  | 0.273465  |       -0.253938 |       -0.169672 | [114032 130880  15335] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |

## digital_twin_properties_bcm_coreg_apl_vess_fwd
schema: digital_twin_properties

description: 
Functional properties computed from the digital twin model from Wang et al. 2025, merged with apl_functional_coreg_vess_fwd.

    cc_abs : Test set performance of the digital twin model unit, higher is better
    cc_max : Neuron variability score used to normalize digital twin model unit performance
    cc_norm : Normalized model unit performance, higher is better
    OSI : orientation selectivity index
    DSI : direction selectivity index
    gOSI : global orientation selectivity index
    gDSI : global direction selectivity index
    pref_ori : Preferred orientation in degrees (0 - 180), vertical bar moving right is 0 and orientation increases counter-clockwise
    pref_dir : Preferred direction in degrees (0 - 360), vertical bar moving right is 0 and orientation increases counter-clockwise
    readout_loc_x : X coordinate of the readout location, an approximation of receptive field center in stimulus space; (-1, -1) top-left, (1, 1) bottom-right
    readout_loc_y : Y coordinate of the readout location, an approximation of receptive field center in stimulus space; (-1, -1) top-left, (1, 1) bottom-right

This table is maintained by Stelios Papadopoulos @ Tolias Lab

For the use of this table please cite:
1) https://doi.org/10.1038/s41586-025-08829-y; Wang et al., 2025
2) https://doi.org/10.1038/s41586-025-08840-3; Ding, Fahey, and Papadopoulos et al., 2025
3) https://doi.org/10.1038/s41586-025-08790-w; MICrONS Consortium et al., 2025
 [Note: This table 'digital_twin_properties_bcm_coreg_apl_vess_fwd' will update the 'target_id' foreign_key when updates are made to the 'nucleus_detection_v0' table] 

columns: ['id_ref', 'created_ref', 'valid_ref', 'volume', 'pt_supervoxel_id', 'pt_root_id', 'id', 'created', 'valid', 'target_id', 'session', 'scan_idx', 'unit_id', 'pref_ori', 'pref_dir', 'gOSI', 'gDSI', 'cc_abs', 'cc_max', 'cc_norm', 'OSI', 'DSI', 'readout_loc_x', 'readout_loc_y', 'pt_position', 'bb_start_position', 'bb_end_position']

|    |   id_ref | created_ref                      | valid_ref   |   volume |   pt_supervoxel_id |         pt_root_id |    id | created                          | valid   |   target_id |   session |   scan_idx |   unit_id |   pref_ori |   pref_dir |      gOSI |      gDSI |   cc_abs |   cc_max |   cc_norm |       OSI |       DSI |   readout_loc_x |   readout_loc_y | pt_position            | bb_start_position   | bb_end_position   |
|---:|---------:|:---------------------------------|:------------|---------:|-------------------:|-------------------:|------:|:---------------------------------|:--------|------------:|----------:|-----------:|----------:|-----------:|-----------:|----------:|----------:|---------:|---------:|----------:|----------:|----------:|----------------:|----------------:|:-----------------------|:--------------------|:------------------|
|  0 |   106251 | 2020-09-28 22:44:30.237033+00:00 | True        |  256.97  |  79946864845341470 | 864691135969663333 | 36042 | 2025-07-30 18:11:18.767569+00:00 | True    |      106251 |         9 |          4 |      2388 |  124.556   |    304.556 | 0.210832  | 0.0323181 | 0.219483 | 0.59398  |  0.369512 | 0.40676   | 0.0263991 |       -0.214859 |       -0.20811  | [110304 120192  15023] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  1 |   106265 | 2020-09-28 22:42:15.032352+00:00 | True        |  217.643 |  80439446054697150 | 864691136130351441 | 36125 | 2025-07-30 18:11:18.825463+00:00 | True    |      106265 |         9 |          4 |      2484 |    2.63202 |    182.632 | 0.0728539 | 0.0595035 | 0.330504 | 0.411639 |  0.802898 | 0.143274  | 0.151337  |       -0.255015 |       -0.163746 | [113952 120208  15112] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  2 |   106431 | 2020-09-28 22:43:14.353734+00:00 | True        |  289.367 |  79947689680746397 | 864691135725228203 | 36169 | 2025-07-30 18:11:18.854936+00:00 | True    |      106431 |         9 |          4 |      2536 |   65.9691  |    245.969 | 0.0231131 | 0.0842573 | 0.4073   | 0.504181 |  0.807845 | 0.0286045 | 0.282779  |       -0.252097 |       -0.143649 | [110288 126400  16890] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  3 |   106436 | 2020-09-28 22:41:20.272045+00:00 | True        |  315.426 |  80018127211293829 | 864691135133547168 | 35902 | 2025-07-30 18:11:18.670003+00:00 | True    |      106436 |         9 |          4 |      2212 |   20.3977  |    200.398 | 0.109068  | 0.0973767 | 0.482162 | 0.719846 |  0.669813 | 0.213756  | 0.276174  |       -0.202574 |       -0.196255 | [110768 126704  17166] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |
|  4 |   106518 | 2020-09-28 22:41:14.258730+00:00 | True        |  309.334 |  80299602255318683 | 864691135210056384 | 36156 | 2025-07-30 18:11:18.846229+00:00 | True    |      106518 |         9 |          4 |      2521 |  103.407   |    283.407 | 0.118369  | 0.0749642 | 0.558106 | 0.830425 |  0.672073 | 0.231545  | 0.257083  |       -0.244328 |       -0.163041 | [112960 126672  17863] | [<NA> <NA> <NA>]    | [<NA> <NA> <NA>]  |

