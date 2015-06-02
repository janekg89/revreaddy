set l [ expr 13. / 2. ]
set s solid
set w 2
set c 16
graphics top color $c
graphics top line "[expr -1. * $l] [expr -1 * $l] [expr -1. * $l]" "[expr -1. * $l] [expr -1. * $l]  $l" width $w style $s
graphics top line "[expr -1. * $l] [expr -1 * $l] [expr -1. * $l]" "[expr -1. * $l]  $l [expr -1. * $l]" width $w style $s
graphics top line "[expr -1. * $l] [expr -1 * $l] [expr -1. * $l]" " $l [expr -1. * $l] [expr -1. * $l]" width $w style $s
graphics top line " $l [expr -1. * $l]  $l" "[expr -1. * $l] [expr -1. * $l]  $l" width $w style $s
graphics top line " $l [expr -1. * $l]  $l" " $l [expr -1. * $l] [expr -1. * $l]" width $w style $s
graphics top line " $l [expr -1. * $l]  $l" " $l  $l  $l" width $w style $s
graphics top line "[expr -1. * $l]  $l  $l" "[expr -1. * $l] [expr -1. * $l]  $l" width $w style $s
graphics top line "[expr -1. * $l]  $l  $l" " $l  $l  $l" width $w style $s
graphics top line "[expr -1. * $l]  $l  $l" "[expr -1. * $l]  $l [expr -1. * $l]" width $w style $s
graphics top line " $l  $l [expr -1. * $l]" "[expr -1. * $l]  $l [expr -1. * $l]" width $w style $s
graphics top line " $l  $l [expr -1. * $l]" " $l  $l  $l" width $w style $s
graphics top line " $l  $l [expr -1. * $l]" " $l [expr -1. * $l] [expr -1. * $l]" width $w style $s
