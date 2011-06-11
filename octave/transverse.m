% Compute transverse flux

% Up flux if signA > 0
% Down flux if signA < 0
function Z = transverseX(signA, lambda, r1, r2, r3, uTilde, vTilde, cTilde)
  if (sign(lambda(1)) == signA)
    A = lambda(1)*alpha(1)*r1;
  end
  if (sign(lambda(2)) == signA)
    A = A + lambda(2)*alpha(2)*r2;
  end
  if (sign(lambda(3)) == signA)
    A = A + lambda(3)*alpha(3)*r3;
  end

  if (sign(lambda(1)) == signA)
    beta1 = ((uTilde+cTilde)*A(1)-A(2))/(2*cTilde);
    Z = beta1*lambda(1)*r1;
  end
  if (sign(lambda(2)) == signA)
    beta1 = -vTilde*A(1)+A(3);
    Z = Z+beta1*lambda(2)*r2;
  end
  if (sign(lambda(3)) == signA)
    beta1 = (-(uTilde-cTilde)*A(1)+A(2))/(2*cTilde);
    Z = Z+beta1*lambda(3)*r3;
  end
end

function Z = transverseY(signA, lambda, r1, r2, r3, uTilde, vTilde, cTilde)
  if (sign(lambda(1)) == signA)
    A = lambda(1)*alpha(1)*r1;
  end
  if (sign(lambda(2)) == signA)
    A = A + lambda(2)*alpha(2)*r2;
  end
  if (sign(lambda(3)) == signA)
    A = A + lambda(3)*alpha(3)*r3;
  end

  if (sign(lambda(1)) == signA)
    beta1 = ((vTilde+cTilde)*A(1)-A(3))/(2*cTilde);
    Z = beta1*lambda(1)*r1;
  end
  if (sign(lambda(2)) == signA)
    beta1 = uTilde*A(1)-A(2);
    Z = Z+beta1*lambda(2)*r2;
  end
  if (sign(lambda(3)) == signA)
    beta1 = (-(vTilde-cTilde)*A(1)+A(3))/(2*cTilde);
    Z = Z+beta1*lambda(3)*r3;
  end
end
