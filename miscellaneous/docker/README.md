The two repositories have docker files that are to be used to generate build and run environments for compass.

Logged into the compass docker regisry:

```shell
docker login registry.gitlab.inria.fr
```

then in the `build` directory:

```shell
docker build -t registry.gitlab.inria.fr/charms/compass/build-environment .
docker push registry.gitlab.inria.fr/charms/compass/build-environment 
```

and in the `run` directory:

```shell
docker build -t registry.gitlab.inria.fr/charms/compass/run-environment .
docker push registry.gitlab.inria.fr/charms/compass/run-environment 
```
