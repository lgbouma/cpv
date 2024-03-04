###############################################################################
Original:  Tue Jun 13 15:44:42 2023
Updates:   Mon Mar  4 10:00:12 2024
Author: luke

This directory contains CPV-relevant literature information.
This mostly means catalogs of previously published CPVs.
There is also a set of Gagne's associations and ages.

These were processed using `drivers/process_literature_cpvs.py`

###############################################################################

------------------------------------------
Stauffer2017 table 1:  K2 Campaign2
Abstract "Total 23 VLM members of ρ Oph and USco; 19 are 'class1' or 'class2',
  other four are RIK-210 analogs."

* N=11 Class 1 = scallop shell
* N=8 Class 2 = persistent short-duration flux dip
* N=4 Class 3 = transient short duration flux dip.  (Includes RIK-210, 5.7 days!), and a few shorter-duration cases...

Class 1 and 2 seem plausible...

Including RIK-210 seems like a stretch, IMO.  Perhaps include?  But w/ a flag.
Else omit...  (probably omit, b/c Stauffer himself omitted these objects in
subsequent papers.  but this is basically inconsistent, b/c identification of
these objects by shape included many "interlopers" in Stauffer2018 and
Stauffer2021 that imho do not meet the criterion of being obvious members.)

-> `get_stauffer17_xmatch` returns 23 objects.

------------------------------------------
Stauffer2018 table 1:
"Class 1" only, in USco, rhoOph, and Taurus

"...We have now identified an additional eight probable members of the same
star-forming region plus three stars in the Taurus star-forming region with
this same light curve morphology"

* N=11 Class 1... except really it's N=12, because EPIC 204060981 has two CPVs.

Yikes.

------------------------------------------
Stauffer2021 table 1:

USco TESS scallops

-----
Stauffer2021 table 2:

USco TESS "class2" persistent flux dips
Includes FWZI's

--------------------
COUNT

Stauffer17 class1+2 = 29
Stuaffer18 class1 += 8
Stauffer21 class1 += 28
Zhan19/Gunther22 += 10
Popinchalk23 += 3
