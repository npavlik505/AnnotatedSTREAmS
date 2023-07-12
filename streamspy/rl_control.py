import math
import random
from collections import namedtuple, deque
from itertools import count
from typing import Deque

# pytorch imports
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F

class RLParameters():
    # batch size: number of transitions we sample from the replay buffer
    # gamma: discount factor
    # epsilon_start: starting value of epsilon
    # epsilon_end: ending value of epsilon (this will be the minimum value of epsilon for long simulations)
    # tau: exponential decay rate factor for epsilon
    # lr: learning rate
    def __init__(self, batch_size: int, gamma: float, epsilon_start: float, epsilon_end: float, epsilon_decay: float, tau: float, lr: float):
        self.batch_size   = batch_size
        self.gamma        = gamma
        self.epsilon_start= epsilon_start
        self.epsilon_end  = epsilon_end
        self.epsilon_decay= epsilon_decay
        self.tau          = tau
        self.lr           = lr

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

        if self.device == "cpu":
            raise ValueError("no GPU available for RL - this should not happen")

    def epsilon_threshold(self, step: int):
        threshold = self.epsilon_end + (self.epsilon_start- self.epsilon_end) * \
            math.exp(-1. * step/ self.epsilon_start)

        return threshold

class Transition():
    def __init__(self, state:float, action: float, next_state: float, reward: float):
        self.state      = state      
        self.action     = action     
        self.next_state = next_state 
        self.reward     = reward     

class ReplayMemory(object):

    def __init__(self, capacity):
        self.memory : Deque = deque([], maxlen=capacity)

    def push(self, *args):
        """Save a transition"""
        self.memory.append(Transition(*args))

    def sample(self, batch_size):
        return random.sample(self.memory, batch_size)

    def __len__(self):
        return len(self.memory)

class DQN(nn.Module):
    def __init__(self, n_observations, n_actions):
        super(DQN, self).__init__()
        self.layer1 = nn.Linear(n_observations, 128)
        self.layer2 = nn.Linear(128, 128)
        self.layer3 = nn.Linear(128, n_actions)

    # Called with either one element to determine next action, or a batch
    # during optimization. Returns tensor([[left0exp,right0exp]...]).
    def forward(self, x):
        x = F.relu(self.layer1(x))
        x = F.relu(self.layer2(x))
        return self.layer3(x)

def select_action(params: RLParameters, step: int, policy_net: DQN, state):
    sample = random.random()

    eps_threshold = params.epsilon_threshold(step)

    if sample > eps_threshold:
        with torch.no_grad():
            # t.max(1) will return the largest column value of each row.
            # second column on max result is index of where max element was
            # found, so we pick action with the larger expected reward.
            return policy_net(state).max(1)[1].view(1, 1)
    else:
        jet_velocity = random.uniform(0., 1.)
        return torch.tensor([[jet_velocity]], device=params.device, dtype=torch.long)
