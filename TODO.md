- Fix stairs behaviour - it currently steps on exactly one element after expected due to the way 'ceil' function rounds up. Fix and uncomment test once done.

- Check PowerLawPlateau is in new format, add if not

- Check if generalised Maxwell model is in new format, add if not

- Update stepgen (with sigmoidal transition), noisegen and repeatdata (with sigmoidal transition) so that they integrate with new workflow

- decide on style guide? frequency_spec vs frequencyspec, time_line vs timeline.

- Remove use of quasinull, replace with Union{Nothing, x} types and remove quasinull function itself

- Keep increasing test coverage

- Add purely FFT based convolution to side-step time domain Mittag-Leffler function

