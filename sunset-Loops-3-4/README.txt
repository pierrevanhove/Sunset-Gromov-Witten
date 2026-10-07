Ancillary files for H. Maennel and P. Vanhove, "The three- and four-loop
equal-mass sunset integrals and local Calabi-Yau fourfolds and fivefolds",
appendix "Reproducible interval certificate".

certify.py, validate.py   the two programs printed in that appendix
                          (byte-identical to the listings)
certificate_32.json, certificate_48.json, validation.json, exact_germ_check.json
                          their outputs, as archived

Produced with Python 3.12.3 and python-flint 0.9.0
(requirements: Python >= 3.10, python-flint==0.9.0).
Regenerate and validate; keep assertions enabled (do not use python3 -O):

  python3 certify.py \
    && python3 certify.py --bits 1500 --order 420 --match 48 \
    && python3 validate.py

A fresh run reproduces validation.json and exact_germ_check.json byte for
byte, and the two certificates in every field except "seconds"
(wall-clock time).

SHA-256 (also in the paper, Table "SHA-256 checksums of the ancillary files"):
f546eab696a7c8dab6e05ef77a54e943658979ccadc512f500d9f110d14ae934  certify.py
9e806a34ca4d6d52f97b0afe44837c94ba98fb502f8295db76b595a53962741d  validate.py
a83d5a9efb323c8828976058ce28bb10f453d88fbc96f905b0d1a153d2c3f7f3  certificate_32.json
ffa6d28a620c73a62e7f221be086c8ab2d2d03b861abeb7b0208b605ae5c9bbb  certificate_48.json
dc849dbe9050a50b980a62be72fa2d4e292c87d525e97e0cd1945de35fd34c0b  validation.json
32c481b01a9741f143408243bd65c772eac7b07342fa2b2d8556de78106247bf  exact_germ_check.json
