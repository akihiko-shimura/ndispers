.. _intro-install:

Installation
============

ndispers works on Python_ 3.7 or higher.

You can install it using pip_::

    $ pip install ndispers

Check the version installed is the latest_. If not, specify the latest version by::

    $ pip install ndispers==x.y.z

To uninstall::

    $ pip uninstall ndispers

To update version::

    $ pip install -U ndispers


Testing your installation
-------------------------

You can check that ndispers is working correctly by starting up a python shell.

.. code-block:: python

    In [1]: import ndispers as nd

    In [2]: nd.media.crystals.BetaBBO_Eimerl1987()
    Out[2]: <ndispers.media.crystals._betaBBO_Eimerl1987.BetaBBO at 0x105368d50>

Here you can see that beta-BBO object is created (Please do not forget that in In[2] :code:`()` is needed to get a mediaum object).
If it fails, feel free to open an issue in `issue tracker`_.

.. note::

    It is strongly reccomended using IPython_ shell or JupyterLab_ because they have `Tab completion`_ function, 
    which is very useful to search for available crystals and glasses.


.. _Python: http://www.python.org/
.. _pip: http://www.pip-installer.org/
.. _latest: https://github.com/akihiko-shimura/ndispers/releases
.. _`issue tracker`: https://github.com/akihiko-shimura/ndispers/issues
.. _IPython: https://ipython.org/documentation.html
.. _JupyterLab: https://jupyter.org/
.. _`tab completion`: https://ipython.readthedocs.io/en/stable/interactive/tutorial.html#tab-completion


