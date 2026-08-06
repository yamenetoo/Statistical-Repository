Tutorial
========

Quick start
-----------

.. code-block:: python

   from unit_transform_dist import make_unit_transform

   dist = make_unit_transform("gamma", base_shape_params=(2.5,))
   dist.pdf(0.3)
   dist.rvs(size=1000, random_state=42)

Fitting and model selection
----------------------------

.. code-block:: python

   from unit_transform_dist import auto_select_distribution

   results = auto_select_distribution(data, criterion="aic", top_n=5)
   best = results[0]["dist"]

Regression on proportions
---------------------------

.. code-block:: python

   from unit_transform_dist import UnitTransformRegression

   model = UnitTransformRegression("gamma", link="log")
   model.fit(X, y)
   model.predict_mean(X_new)

See the README for the full feature list, including mixtures, Bayesian
estimation, censored/truncated fitting, and diagnostics.
