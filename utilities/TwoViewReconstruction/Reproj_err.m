function err = Reproj_err(pix_meas, p)

pix_est = p/p(3);
err = norm(pix_meas - pix_est);