
# Où se trouve mpi ?

* Case
  - `_prepare_distribute`  -> gather ?
  - `_distribute_mesh`  -> scatter
  - `_distribute_properties`  -> scatter
  - `_distribute_zones`  -> scatter
* ZoneManager
  - `BaseManager._new_id`  -> bcast
  - `BaseManager._id_from_parts`  -> gather, bcast
  - `Calculator._bool_all`  -> allreduce (all)
  - `Calculator._bool_any`  -> allreduce (any)
* Property
  - rien ?


# Idée API Case

Rien de définitif pour l'instant mais voici une ébauche:

```python
case = Case.from_mesh_spec(...)
case.mesh

case.positions.cells -> ArrayProperty (x, y, z)
case.location.cells -> ZoneManager sur cells
case.prop1

mapping_proc = overlap_strategy(metis(np.mesh))

case.distribute(mapping_proc)
case.is_distributed

result = ComPASS.solve(case)
ComPASS.solve(case)  # -> case.solution
```


# Zone API

* on pourrait préférer `zone.where()` à `zone.content`
  - `()` pour marquer le fait qu'un calcul est nécessaire
  - `where()` en référence à `np.where(arr)` qui renvoit un array d'indices

# Où et quand indique-t-on la dimension de l'espace ?

* pour l'instant, on peut faire Property(Tensor(2)) mais il faut indiquer la dim (2) à la création.
  - c'est dommage car on souhaite aussi mettre ces properties dans case et celui doit pouvoir s'adapter à la dimension du problème.
* et dans case, ça va où ?
