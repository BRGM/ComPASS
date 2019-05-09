
z = objet_qui_connait_tout.zones()
v = objet_qui_connait_tout.variables()


# pour mettre une valeur partout
# (choisir une syntaxe)
v.porosity = 0.23
v.porosity[...] = 0.23
v.porosity.everywhere = 0.23


centre = [2.3e3, 3.4e3, -100]
rayon = 25

xyz = v.locations[z.cells]  # positions des centres des cellules
sable = z.cells.from_mask(
    np.linalg.norm(xyz - centre) < rayon
)

v.porosity[sable] = 0.4


# conversion de types compliqués
v.permeability[z.cells] = 2  # tenseur isotrope
v.permeability[z.cells] = [[3, 0], [4, 0]]


## récupération des arrays
##

# array porté par z.cells.content()
v.porosity[z.cells]


###################################################
###################################################


## sur un maillage déjà colorié

"""
z.cells.manager.get(1)
renvoie la zone des cellules de tag valant 1 dans le maillage
si le tag n'existe pas, renvoie la zone vide

v.porosity[z1, z2] est équivalent à
  v.porosity[z1]
  v.porosity[z2]
si z1 et z2 sont dans le même manager, on peut aussi faire
  v.porosity[z1 | z2]
mais pas si z1 et z2 ont des managers différents.
syntaxe pour réduire des répétitions

"""

z = objet_qui_connait_tout.zones()
v = objet_qui_connait_tout.variables()

cells = z.cells.manager.get
nodes = z.nodes.manager.get
faces = z.faces.manager.get

v.porosity[cells(1), nodes(1)] = 0.23
v.permeability[cells(1)] = 5e-3

v.porosity[cells(2), nodes(2)] = 0.23
v.permeability[cells(2)] = [[3e-3, 0], [0, 1e-3]]


###############################
###############################

""" objet_qui_connait_tout.distribute()

Zict est overkill: Property n'utilise que __getitem__/__setitem__.
Zone est hashable pour pouvoir aller dans Zict.
Zone hashable est en contradiction avec distribute().
=> exit Zict: ZoneMap

ZoneMap(zone_manager)
ZoneMap.get(zone)  -> Value or None
ZoneMap.set(zone, value)


0 - le zone manager de chaque proc reçoit les indices de son nouveau domaine
    (own+ghost)

1 - accorder le partitionnement du zone manager en accord avec son nouveau
    domaine: il suffit de construire la zone correspondante

2 - 


"""

###############################
###############################

""" Manitou

Manitou -> Case

"""

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


###############################
###############################

""" Distribution

distribution des zones passe par le manager
 1/ master calcul mng._partition pour chaque proc
 2/ master calcul les zone._block_ids pour chaque proc
 3/ synchro

distribution des property
  * property._data[mng] -> {zone: val, ...}
 1/ calcul de val_new qui remplacera val
 2/ distribution de mng
 3/ maj de property._data[mng] avec les val_new



"""
