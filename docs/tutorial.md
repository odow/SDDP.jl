# Tutorial

Moa.jl is a library for solving Multi-stage stochastic optimisation problems. A classic example in the literature is the so-called "newsvendor problem." However, for this tutorial, we shall adapt it to the "ice-cream salesperson problem."

In this problem, a seller of ice-creams must decide the quantity of ice-cream to purchase over the course of one week during summer.

The seller purchases a quantity of ice-cream at the start of each day from a supplier, and then spends the rest of the day travelling around selling ice-creams. If there is ice-cream that is unsold at the end of the day, it can be stored overnight in a freezer, but incurs a cost of storage.

Furthermore, the price at which the seller purchases the bulk ice-cream varies over time during the week as the supplier seeks to maximise their own profit by charging more during the weekends.

In addition, if the seller knew the demand for ice-cream ahead of time, they could perfectly optimise their purchases to minimise the total cost (or maximise the total profit). However, demand for ice-cream varies stochastically depending on the weather. Therefore, the task of deciding the quantity of ice-cream to purchase is not straight forward and so the seller decides to use Moa.jl to optimise their purchasing decisions.
