import numpy as np
import timeit
import torch
import torch.nn as nn
import torch.optim as optim
from torchdiffeq import odeint_adjoint as odeint

true_y0 = torch.tensor([[2., 0.]])
t = torch.linspace(0., 1.5, 30)
true_A = torch.tensor([[-0.1, 2.0], [-2.0, -0.1]])

class Lambda(nn.Module):
    def forward(self, t, y):
        return torch.mm(y**3, true_A)

with torch.no_grad():
    true_y = odeint(Lambda(), true_y0, t, method='dopri5')

class ODEFunc(nn.Module):

    def __init__(self):
        super(ODEFunc, self).__init__()

        self.net = nn.Sequential(
            nn.Linear(2, 50),
            nn.Tanh(),
            nn.Linear(50, 2),
        )

        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                nn.init.constant_(m.bias, val=0)

    def forward(self, t, y):
        return self.net(y**3)

func = ODEFunc()
optimizer = optim.Adam(func.parameters(), lr=0.05)

def time_func():
    for itr in range(1, 501):
        optimizer.zero_grad()
        pred_y = odeint(func, true_y0, t)
        loss = torch.sum((pred_y - true_y)**2)
        loss.backward()
        optimizer.step()
        with torch.no_grad():
                pred_y = odeint(func, true_y0, t)
                print(torch.sum((pred_y - true_y)**2))


time_func()

func = ODEFunc()
optimizer = optim.Adam(func.parameters(), lr=0.05)
def time_func():
    for itr in range(1, 501):
        optimizer.zero_grad()
        pred_y = odeint(func, true_y0, t)
        loss = torch.sum((pred_y - true_y)**2)
        loss.backward()
        optimizer.step()

timeit.Timer(time_func).timeit(number=1) # 288.965871299999 seconds

with torch.no_grad():
    pred_y = odeint(func, true_y0, t)
    print(torch.sum((pred_y - true_y)**2))

# tensor(0.0596)