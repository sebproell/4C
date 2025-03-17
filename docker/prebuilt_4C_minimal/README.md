>**Disclaimer:** this docker image is experimental. There will be modifications in the future.

This `Dockerfile` creates an docker image with a ready to use 4C executable. Compared to
`prebuilt_4C`, this image only contains the bare minimum to run 4C, so, no source, build folder,
tests, or other testing infrastructure.

This image is for you if you just want a 4C executable or want to use it as image in GitHub Actions.

The 4C executable is in `/home/user/4C/bin/4C` but can also be directly called with `4C` from
anywhere via a symlink.
