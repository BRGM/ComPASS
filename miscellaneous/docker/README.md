These repositories have docker files that are used to generate 
several environments to work with ComPASS:
  - *build-environment* is designed to compile ComPASS,
  - *run-environment* is designed to run ComPASS, the compiled version of the
ComPASS modules are generated using the *build-environment* and installed into
speicific version of this container which are designed to be as slim as possible, 
these images are the ones that are generated at each commit and can be pulled
from the registry (cf. `Dockerfile` in the root directory of ComPASS source),
  - *work-environment* is designed to compile, modify and run ComPASS from a
local file system that has to be mounted on the `localfs` volume of the container.
This container also ships `scipy` and `matplotlib` to run and exploit simulations
efficiently.

To manage permission issues we use the [trick from Deni Bertovic](https://denibertovic.com/posts/handling-permissions-with-docker-volumes/)
 that is implemented in the *switch-user* container.

## Generation of the images

Log into the ComPASS docker regisry:

```shell
docker login registry.gitlab.inria.fr
```

then:

```shell
docker build -t registry.gitlab.inria.fr/charms/compass/base-environment base
docker push registry.gitlab.inria.fr/charms/compass/base-environment 
docker build -t registry.gitlab.inria.fr/charms/compass/build-environment build
docker push registry.gitlab.inria.fr/charms/compass/build-environment 
docker build -t registry.gitlab.inria.fr/charms/compass/run-environment run 
docker push registry.gitlab.inria.fr/charms/compass/run-environment 
docker build -t registry.gitlab.inria.fr/charms/compass/work-environment work
docker push registry.gitlab.inria.fr/charms/compass/work-environment 
```

