# Ole Sigmund's 99 line topology optimization script ported to Julia
function top(nelx, nely, volfrac, penal, rmin)
  # INITIALIZE
  x = fill(volfrac, (nely,nelx))
  dc = zeros(nely, nelx)
  loop = 0
  change = 1.0

  # START ITERATION
  while change > 0.01
    loop += 1
    xold = x

    # FE-ANALYSIS
    U = FE(nelx, nely, x, penal)

    # OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    KE = lk()
    c = 0.0
    for ely = 1:nely
      for elx = 1:nelx
        n1 = (nely + 1) * (elx - 1) + ely
        n2 = (nely + 1) * elx + ely
        edof = [
          2n1 - 1;
          2n1;
          2n2 - 1;
          2n2;
          2n2 + 1;
          2n2 + 2;
          2n1 + 1;
          2n1 + 2;
        ]
        Ue = U[edof]
        compliance = (Ue' * KE * Ue)[1]
        c = c + x[ely,elx]^penal * compliance
        dc[ely,elx] = -penal * x[ely,elx]^(penal - 1) * compliance
      end
    end

    # FILTERING OF SENSITIVITIES
    dc = check(nelx, nely, rmin, x, dc)

    # DESIGN UPDATE BY THE OPTIMALITY CRITERIA METHOD
    x = OC(nelx,nely,x,volfrac,dc)

    # PRINT RESULTS
    change = maximum(abs(x - xold))
    println(@sprintf("It.:  %4i", loop))
    println(@sprintf("Obj.: %10.4f", c))
    println(@sprintf("Vol.: %6.3f", sum(sum(x)) / (nelx*nely)))
    println(@sprintf("ch.:  %6.3f", change))

    # PLOT DENSITIES
    # colormap(gray); imagesc(-x); axis equal; axis tight; axis off;pause(1e-6)
  end

  showall(x)
end

# OPTIMALITY CRITERIA UPDATE
function OC(nelx, nely, x, volfrac, dc)
  l1, l2, move = 0, 100000, 0.2
  xnew = zeros(nelx, nely)

  while (l2 - l1 > 1e-4)
    lmid = (l2 + l1) / 2
    xnew = max(0.001, max(x - move, min(1.0, min(x + move, x .* sqrt(-dc ./ lmid)))))
    if sum(sum(xnew)) - volfrac * nelx * nely > 0
      l1 = lmid
    else
      l2 = lmid
    end
  end

  return xnew
end

# MESH-INDEPENDENCY FILTER
function check(nelx, nely, rmin, x, dc)
  dcn = zeros(nely, nelx)
  radius = floor(Integer, rmin)

  for i in 1:nelx
    for j in 1:nely
      sum = 0.0

      for k in max(i - radius, 1):min(i + radius, nelx)
        for l in max(j - radius, 1):min(j + radius, nely)
          fac = rmin - sqrt((i - k)^2 + (j - l)^2)
          sum = sum + max(0, fac)
          dcn[j,i] = dcn[j,i] + max(0, fac) * x[l,k] * dc[l,k]
        end
      end

      dcn[j,i] = dcn[j,i] / (x[j,i] * sum)
    end
  end

  return dcn
end

# FE-ANALYSIS
function FE(nelx, nely, x, penal)
  KE = lk()
  n = 2 * (nelx + 1) * (nely + 1)
  K = spzeros(n, n)
  F = spzeros(n, 1)
  U = zeros(n, 1)

  for elx in 1:nelx
    for ely in 1:nely
      n1 = (nely + 1) * (elx - 1) + ely
      n2 = (nely + 1) * elx + ely
      edof = [
        2n1 - 1;
        2n1;
        2n2 - 1;
        2n2;
        2n2 + 1;
        2n2 + 2;
        2n1 + 1;
        2n1 + 2;
      ]
      K[edof,edof] = K[edof,edof] + x[ely,elx]^penal * KE
    end
  end

  # DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
  F[2,1] = -1
  fixed_dofs = vcat(1:2:2*(nely + 1), [n])
  all_dofs = collect(1:n)
  free_dofs = setdiff(all_dofs, fixed_dofs)

  K[free_dofs,free_dofs] \ F[free_dofs,:]

  # SOLVING
  U[free_dofs,:] = K[free_dofs,free_dofs] \ F[free_dofs,:]

  return U
end

# ELEMENT STIFFNESS MATRIX
function lk()
  E = 1.0
  nu = 0.3
  k = E / (1 - nu^2) * [
    1/2 - nu/6
    1/8 + nu/8
    -1/4 - nu/12
    -1/8 + 3nu/8
    -1/4 + nu/12
    -1/8 - nu/8
    nu/6
    1/8 - 3nu/8
  ]

  KE = [
    k[1] k[2] k[3] k[4] k[5] k[6] k[7] k[8];
    k[2] k[1] k[8] k[7] k[6] k[5] k[4] k[3];
    k[3] k[8] k[1] k[6] k[7] k[4] k[5] k[2];
    k[4] k[7] k[6] k[1] k[8] k[3] k[2] k[5];
    k[5] k[6] k[7] k[8] k[1] k[2] k[3] k[4];
    k[6] k[5] k[4] k[3] k[2] k[1] k[8] k[7];
    k[7] k[4] k[5] k[2] k[3] k[8] k[1] k[6];
    k[8] k[3] k[2] k[5] k[4] k[7] k[6] k[1];
  ]

  return KE
end
