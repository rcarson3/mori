import numpy as np
import math

quat = np.zeros(4)
vec = np.zeros(3)
rvec = np.zeros(3)

quat[0] = 0.5
quat[1:4] = -0.5
print(quat)

vec[2] = 1.0  

q02 = 2.0 * quat[0]
#//(q_0^2 - ||q||^2)
tmp = (quat[1] * quat[1]) + (quat[2] * quat[2]) + (quat[3] * quat[3])
print(tmp)
q02_m_nq = quat[0] * quat[0] - tmp
print(q02_m_nq)
dot_prod2 = 2.0 * (quat[1] * vec[0] + quat[2] * vec[1] + quat[3] * vec[2])
cross_prod = np.zeros(3)

cross_prod[0] = -quat[3] * vec[1] + quat[2] * vec[2]
cross_prod[1] = quat[3] * vec[0] - quat[1] * vec[2]
cross_prod[2] = -quat[2] * vec[0] + quat[1] * vec[1]

print(cross_prod)
print(dot_prod2)

rvec[0] = vec[0] * q02_m_nq + cross_prod[0] * q02 + quat[1] * dot_prod2
rvec[1] = vec[1] * q02_m_nq + cross_prod[1] * q02 + quat[2] * dot_prod2
rvec[2] = vec[2] * q02_m_nq + cross_prod[2] * q02 + quat[3] * dot_prod2 

print(rvec)