- Check stairs behaviour is as expected - it currently steps on exactly one element after one might expect due to the way the 'ceil' function rounds up. Decide, fix if necessary and and uncomment test once done.

- In definitions.jl, consider replacing empty_rheodata_vector with simply RheoFloat[] as it is self explanatory and actually shorter, more readable. What is the advantage of empty_rheodata_vector? If it must be kept then define it as a 'const' in the same way that RheoFloat is now defined.

- Check PowerLawPlateau is in new format, add if not

- Check if generalised Maxwell model is in new format, add if not

- - Update stepgen (with sigmoidal transition), noisegen and repeatdata (with sigmoidal transition) so that they integrate with new workflow

- Remove use of quasinull, replace with Union{Nothing, x} types and remove quasinull function itself

- Keep increasing test coverage

- Add purely FFT based convolution to side-step time domain Mittag-Leffler function

