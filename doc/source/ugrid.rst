ugrid
++++++++

.. automodule:: ugrid

ugrid (Unstructured grid) is an abstract class for representing the geometry information of the unstructured grid. There are mainly three classes in ugrid module, that is `ugrid_set`, `ugrid_map` and `ugrid_dat`. 

ugrid_set
=========

`ugrid_set` defines the set of unstructured grid, such as the variables in the nodes or the cells.
Its componments and construction method are defined as follows.

.. autoclass:: ugrid_set
   :show-inheritance:
   :members:

To define sets for the specific unstructured grid, the user could use following code:

ugrid_map
=========

`ugrid_map` defines the map between 

.. autoclass:: ugrid_map
   :show-inheritance:
   :members:

ugrid_dat
=========

`ugrid_dat` defines variables associated with set. 

.. autoclass:: ugrid_dat
   :show-inheritance:
   :members: