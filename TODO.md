- Fix 'stairs' loading generator and uncomment test. (It currently steps on exactly one element after it should due to the way the 'ceil' function rounds up).

- Check PowerLawPlateau is in new format, add if not

- Check if generalised Maxwell model is in new format, add if not

- Remove use of quasinull, replace with Union{Nothing, x} types and remove quasinull function itself

- Keep increasing test coverage

- Add purely FFT based convolution to side-step time domain Mittag-Leffler function

- Update stepgen (with sigmoidal transition), noisegen and repeatdata (with sigmoidal transition) so that they integrate with new workflow
