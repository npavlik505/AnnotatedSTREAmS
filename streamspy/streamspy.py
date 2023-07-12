from gymnasium import Env
from typing import Tuple, Any, Optional

ActType = Any
ObsType = Any
RenderFrame = Any

class StreamsSolver(Env):
    def __init__(self):
        pass

    def step(self, action: ActType) -> Tuple[ObsType, float, bool, bool, dict]:
        #print(action)
        pass

    def reset(self, *, seed: int | None = None, options: dict[str, Any] | None = None) -> Tuple[ObsType, dict[str, Any]]:
        pass

    def render(self) -> RenderFrame | list[RenderFrame] | None:
        return None

    def close(self) -> None:
        # solver closes naturally, might close it with MPI
        return None
