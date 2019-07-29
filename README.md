# d18O-ivc-sw
Compute oxygen-isotope chronologies for seawater from foraminifera in sediment cores: ice-volume and temperature corrected

This library will help you go through the following steps:
* Start with raw &delta;<sup>18</sup>O (VPDB) values
* Convert them to VSMOW
* Using Mg/Ca ratios, compute paleotemperatures
* Interpolate missing paleotemperatures
* Using these paleotemperatures, convert the &delta;<sup>18</sup>O (VSMOW) values into &delta;<sup>18</sup>O_<sub>sw</sub>. This subscript ("sw") stands for "seawater", as this step corrects for temperature effects on the isotope ratio within the formainiferal shells and returns the expected seawater &delta;<sup>18</sup>O value.
* Using a sea-level curve, this corrects for ice-volume effects on the bulk ocean &delta;<sup>18</sup>O, thus providing an ice-volume-corrected &delta;<sup>18</sup>O_<sub>ivc-sw</sub>
Error is propagated through each of these steps.

The resultant &delta;<sup>18</sup>O_<sub>ivc-sw</sub> can be used in isotope-mixing studies for paleohydrology.

*Note to users: there is still a bit hard-coded in here; I will fix that eventually but let me know if you need it sooner.*

This process is based on the steps taken by:

**Wickert, A. D., J. X. Mitrovica, C. Williams, and R. S. Anderson (2013), Gradual demise of a thin southern Laurentide ice sheet recorded by Mississippi drainage, *Nature*, 502(*7473*), 668–671, doi:10.1038/nature12609.**
