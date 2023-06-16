import numpy as np
periods = [
    0.7732752056726647,
    0.7731461580225912,
    0.7732813973828846,
    0.7735711202923402,
    0.7729424491300968,
    0.7731911738307666,
    0.7727680841024896,
]
mean_period = np.mean(periods)

print(f"<P> = {mean_period*24:.5f} hr")
print(f"max(P) = {max(periods)*24:.5f} hr")
print(f"min(P) = {min(periods)*24:.5f} hr")
