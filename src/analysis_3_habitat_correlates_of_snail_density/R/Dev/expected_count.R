
stems <- 0
veg <- 0
lake <- 1


linpred_pq <- 0.75 + -0.14 * stems + -0.88 * veg + 2.15 * lake
linpred_cq <- -2.56 + 0.06 * stems + 1.27 * veg - 2.02 * lake

zero_prob_q = 1/ (1 + exp( - linpred_pq ) )

exp(linpred_cq * zero_prob_q) - 0.1757492
