H = [
[1,1,0,1,1,1,1,1,1,0,0,1,0,0,1,0,1,0,0,1,1,0,1,1,1,0,1,1,1,0,0,0,1,0,1,1,0,1,0,0,0,0,0,1,0,0,1,1,0,1,1,0,1,0,1,1,0,1,1,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,0,1,0,1,1,0,0,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,1,0,0,1,0,0,1],
[0,1,1,1,0,0,0,0,0,1,1,0,1,0,0,1,1,1,1,1,1,0,0,1,0,1,0,1,0,0,1,1,1,1,0,1,1,1,0,0,0,1,0,0,0,1,1,1,0,1,1,0,0,1,0,1,0,1,1,0,0,1,1,1,1,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,1,1,0,1,0,1,0,0,1,0,0,0,0,1,1,1,0,1,1,1,1,0,0,1,1,0,0,0,1,1,1,1,0,1,0,0,1,0,0,1,1,1,0,1,0],
[1,0,0,1,0,0,1,1,1,1,1,0,0,1,1,1,1,0,1,0,1,1,0,1,1,0,0,1,0,0,0,0,0,0,1,0,0,1,1,1,1,0,1,0,1,0,1,0,1,0,1,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,1,0,1,1,1,1,1,1,0,0,0,1,1,1,0,0,1,0,1,1,1,0,0,0,0,0,0,0,1,0,1,1,0,0,1,1,1,0,0,0,1,0,1,1,1,1,0,0,0,0,1,1,1,0,1,1,1,1,1,1,0,0,0,0],
[1,1,0,0,1,1,0,1,0,0,1,1,0,1,1,1,0,0,1,0,0,0,0,1,0,1,0,1,1,0,0,1,1,0,1,0,1,1,0,0,0,1,1,1,0,1,1,0,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,0,1,0,0,1,0,0,1,1,1,0,1,1,1,0,1,0,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,1,1,0,0,0,1,1,1,1,1,1,0,0,0,1,0,1,0,0,0,1,0,0,0,1,1,1,0,0,0],
[1,0,0,0,1,1,1,0,1,1,1,0,0,1,1,0,1,1,0,1,1,1,0,1,1,1,0,0,0,1,0,0,1,1,0,1,1,0,0,1,1,0,1,1,0,1,0,1,0,1,1,0,0,1,0,1,0,1,1,0,1,0,1,1,0,0,1,0,0,1,1,1,1,0,1,1,1,0,0,0,0,1,1,1,0,0,0,0,1,1,1,0,0,0,1,1,1,1,1,1,0,0,0,0,1,1,1,0,1,1,1,1,0,1,1,1,1,0,1,1,0,1,0,0,0,1,1,1,0,0,0,1,0,1,0],
[0,0,1,1,0,0,1,0,0,1,1,1,1,0,1,1,1,0,1,0,1,0,0,0,1,1,1,1,1,1,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,1,1,1,1,0,1,1,1,1,1,0,1,1,0,0,0,0,1,0,1,0,0,0,1,0,1,0,1,0,1,0,0,1,1,1,1,1,0,0,0,1,0,0,1,0,0,1,1,1,1,1,1,1,0,0,1,0,0,1,1,1,0,1,1,1,1,0,0,0,0,0,0,1,0,1,0,1,1,0,0,0,1,0,0,1,1,1,1,0],
[1,0,1,1,1,1,1,0,0,0,0,0,1,0,0,1,0,1,0,0,1,1,0,1,0,1,0,0,0,1,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,0,0,0,1,0,0,0,0,1,1,0,1,0,0,1,1,1,1,0,0,0,1,1,1,0,1,0,1,0,1,1,1,0,0,0,1,0,1,1,0,1,1,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,1,1,0,1,1,0,1,1,1,0,1,1,1,0,0,1,1,0,0,0,0,0,1,0,0,0],
[0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,0,1,0,1,0,1,1,1,1,1,0,1,1,1,0,1,0,1,1,0,1,1,1,0,0,1,0,1,1,0,1,1,1,0,0,1,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,1,0,1,0,1,1,0,0,1,1,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,0,1,0,0,0,1,1,1,0,1,0,0,1,0,0,1,1,0,0,0,1,0,1,0,1,0,1,0,0,0,1,0,1,1,1],
[1,1,1,1,0,1,0,1,0,1,0,1,1,1,1,0,1,1,0,1,1,0,1,1,0,0,0,1,1,0,0,0,0,1,0,1,0,1,1,0,1,1,0,0,1,1,0,0,0,0,1,0,0,1,0,1,0,0,0,1,0,1,1,1,1,1,0,1,0,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,0,0,1,1,1,1,1,1,0,0,1,0,0,0,0,1,0,0,1,0,1,1,0,0,1,0,1,0,1,1,1,1,1,1,0,0,1,1,0,1,0,1,0,0,0,1,0,0,1],
[0,1,1,1,1,0,0,1,1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,0,0,0,0,0,0,1,0,1,1,0,0,0,1,1,0,0,1,0,0,1,1,1,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,1,1,0,0,1,0,1,0,0,1,1,0,1,1,0,1,1,0,1,0,1,1,0,1,1,1,0,1,0,1,0,1,1,0,0,0,1,0,1,1,1,0,1,1,0,0,1,1,1,0,1,0,1,1,0,1,1,0,0,1,0,1,1,1,0,1,1,0,1],
[1,1,0,1,0,0,0,1,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,0,1,1,0,1,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,1,1,1,0,1,1,1,1,0,0,0,0,0,0,1,0,1,1,1,0,0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,1,1,1,1,0,0,1,0,1,0,1,0,0,1,1,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,1,1,0,1,1,1,0,0,0,0],
[1,1,1,0,0,1,1,1,0,0,1,1,0,1,1,0,0,1,0,1,1,1,1,0,0,1,1,1,0,1,1,1,1,1,1,1,0,1,1,0,1,0,1,0,0,1,1,0,0,1,1,1,1,0,1,1,1,0,1,0,0,1,1,1,0,1,0,0,1,0,1,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,1,1,1,0,1,0,1,0,1,0,1,1,1,1,0,0,1,0,0,1,0,0,0,1,0,1,0,1,0,0,1,0,1,1,0,1,1,0,1,1],
[1,0,1,0,1,0,1,1,1,0,1,0,1,1,0,1,0,1,1,1,0,0,1,1,1,1,0,0,0,1,0,0,0,1,0,1,0,1,1,0,1,1,0,0,1,0,1,0,1,0,1,0,1,0,1,0,0,1,1,0,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,1,1,1,1,1,1,0,0,0,0,0,1,1,0,1,0,0,1,0,0,1,1,0,0,0,1,1,1,0,0,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,0,1,1,1,0,1,0,0,0,1,1,0,1],
[1,1,1,0,0,1,0,0,0,0,1,0,0,1,1,1,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,0,1,0,1,0,0,0,1,0,0,1,0,0,0,0,1,0,1,1,0,1,1,0,0,1,1,0,1,0,0,1,1,1,0,0,0,1,0,0,1,0,1,1,0,1,1,0,0,1,1,0,1,0,1,1,0,0,0,1,0,1,1,1,0,0,0,1,0,1,1,1,0,1,1,1,1,0,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,0,1,1,0,0,1,1,1,1,1,1,0],
[0,1,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,1,1,0,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,1,0,1,0,0,1,1,1,1,1,0,0,1,0,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,0,1,0,1,0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,1,1,1,1,0,1,0,0,0,0,1,1,0,1,0,1,0,0,1],
[0,1,0,0,0,1,0,0,0,1,1,1,0,0,1,0,1,1,1,1,1,1,0,1,0,0,1,1,1,1,0,1,0,0,0,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,1,1,1,0,1,0,0,1,0,0,0,1,0,0,1,0,1,1,1,0,1,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,1,0,0,1,0,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,1,1,0,0,0,0,0,1,0,0,1,1,1,1,1,1,1,1,0],
[0,1,1,1,0,0,0,1,0,1,1,0,0,0,0,1,0,1,0,1,1,1,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,1,0,1,1,1,1,0,1,1,0,1,1,1,0,0,1,0,1,1,0,1,0,1,0,0,1,0,0,1,1,0,0,1,0,0,1,0,1,1,1,0,1,1,1,0,0,0,0,0,1,0,1,1,1,1,1,1,0,0,0,1,0,1,1,1,1,1,1,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,1,1,1,0,1,0],
[0,0,1,0,0,0,1,1,0,1,0,1,0,0,0,0,0,1,0,0,1,0,1,0,0,1,1,1,0,1,0,0,1,1,1,1,1,0,1,1,0,0,1,0,0,0,1,0,1,1,0,0,0,1,1,1,1,0,1,1,1,0,0,0,1,1,1,0,1,0,1,0,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,1,1,1,1,0,1,1,0,0,0,1,0,0,0,0,0,1,1,0,1,1,1,1,1,0,0,1,1,0,0,0,1,0,0,1,1,0,1,0,0,1],
[1,0,1,0,0,1,1,1,1,0,0,1,0,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,1,0,1,0,0,1,0,0,0,0,1,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,0,0,0,1,0,1,1,1,1,1,0,1,0,0,0,0,0,0,1,1,0,1,1,0,1,1,0,1,0,1,1,0,1,1,0,1,0,1,0,1,1,0,1,1,1,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,0,0,1,1,0,0,1,1,1],
[0,1,1,0,0,0,1,1,1,0,1,1,1,1,1,0,0,1,0,0,0,1,1,0,1,0,1,0,1,0,1,1,0,1,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,0,0,0,0,1,0,1,0,1,0,0,1,0,0,0,0,1,0,0,1,1,0,1,0,1,0,0,0,1,0,1,0,1,1,0,1,1,0,1,0,1,0,1,1,1,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,0,0,1,0,0,1,1,0,0,1,1,1,1,1,0,0,0,0,1,0,0,1,1,0,0,1],
[0,1,1,0,0,1,1,1,1,0,0,0,1,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,1,0,0,1,0,1,0,1,0,0,1,0,0,0,0,0,0,0,1,1,0,1,1,1,1,0,0,1,0,0,0,1,1,1,1,1,0,1,1,1,0,0,1,0,1,1,0,1,0,1,0,0,0,1,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,1,1,0],
[0,0,1,1,1,1,1,0,0,0,0,1,0,0,0,0,1,1,0,1,1,0,0,0,1,1,0,0,1,0,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,1,1,0,1,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,0,1,1,0,1,0,0,0,1,0,1,0,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,1,0,0,1,1,1,0,0,1,1,1,0,1,0,0,0,1,0,0,0,1,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,1,0,1,1,1,0],
[0,1,0,1,0,0,1,1,1,1,1,1,0,0,1,0,1,1,0,0,0,0,1,0,1,1,1,0,1,1,1,1,0,0,1,1,1,0,1,0,0,0,0,1,0,0,0,1,1,0,0,0,1,0,1,0,1,1,0,1,1,0,1,0,0,0,1,1,1,0,0,1,1,1,0,0,1,1,1,1,1,1,0,1,0,1,0,1,0,0,1,1,1,1,1,1,0,0,1,0,1,0,1,0,0,1,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,1,0,1],
[1,1,0,0,0,0,0,0,1,0,1,0,0,0,1,1,0,1,1,1,0,1,1,0,0,0,0,1,1,0,1,1,1,1,0,0,1,0,1,0,0,0,0,1,1,0,1,0,0,1,0,0,0,1,1,1,1,1,1,0,0,1,0,0,1,1,0,0,0,1,1,1,1,1,0,1,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,1,0,1,0,0,0,1,0,1,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,1,1,1,1,0,0,1],
[0,0,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,1,1,1,0,1,1,0,0,1,0,0,0,1,1,0,1,1,0,1,0,1,1,1,0,1,0,1,1,0,0,1,1,0,0,0,1,1,0,1,1,1,1,0,1,0,1,0,0,1,1,1,0,1,1,0,1,0,1,0,1,0,0,0,1,1,0,1,1,0,1,0,0,0,0,0,0,1,0,0,1,1,1,1,1,1,0,1,1,1,0,1,0,1,1,0],
[1,1,0,0,1,0,0,1,1,0,1,1,0,1,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,1,1,0,1,0,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,0,1,1,1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,1,1,0,0,1,1,0,0,0,0,0,0,1,1,1,1,0,0,1,1,0,1,0,1,0,1,0,1,1,0,1,1,0,0,0,1,1],
[0,0,0,0,0,1,0,1,0,0,0,1,1,0,1,0,0,0,1,1,1,0,0,1,0,1,0,1,1,1,1,0,1,0,1,1,0,1,0,0,0,0,1,0,1,0,1,1,0,0,1,0,1,0,1,0,0,1,1,0,0,0,1,1,0,1,1,0,0,0,1,0,1,1,1,1,0,0,1,1,0,0,1,1,1,1,1,1,0,1,0,1,0,1,0,0,0,1,0,1,0,1,1,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1],
[0,0,0,1,1,1,0,1,1,0,1,0,0,0,1,1,1,1,1,0,1,1,0,1,1,1,1,0,0,1,1,1,0,0,1,0,1,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,1,1,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0,1,1,1,0,1,0,1,1,1,0,1,1,1,1,1,0,1,1,0,1,1,0,0,1,0,0,0,0,0,1,1,0,0,1,0,1,1,1,0,1,1,1,1],
[1,0,1,0,1,0,0,1,1,1,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,1,1,0,1,1,0,1,1,1,1,0,1,0,1,0,1,1,1,1,0,1,0,0,1,0,0,0,1,1,1,0,1,0,0,1,1,0,0,1,0,1,0,1,0,0,0,0,1,1,1,0,1,1,1,1,1,0,0,0,1,0,1,1,0,0,0,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,1,1,1,0,0,0],
[1,1,1,0,1,1,0,1,0,1,1,0,0,0,0,0,1,0,0,1,1,1,1,1,0,1,1,1,0,1,0,1,0,1,1,0,1,1,1,1,0,1,0,0,0,1,1,0,0,1,0,0,1,1,0,0,1,1,0,1,0,0,1,1,1,1,1,1,0,1,0,0,1,0,1,0,1,0,0,1,0,1,1,1,1,0,1,0,1,0,1,0,1,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,0,0,1,1,1,1,0,1,1,1,1,0,1,0,0,1,1,0,1,0,1,1,0,1,1,0,0],
[0,1,0,0,1,0,0,0,0,0,1,1,1,0,0,1,1,0,1,1,1,1,0,0,0,1,1,1,1,0,0,1,1,0,1,1,0,0,0,0,1,1,0,1,1,1,0,0,0,0,1,1,0,0,0,0,1,0,1,0,1,0,1,0,0,0,1,0,0,0,1,1,0,0,0,1,1,0,0,1,1,1,0,1,0,0,0,1,1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,0,1,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,1,0,1,1,0,1,1,1,0,1,0,0,1,1,1],
[1,1,0,0,1,0,0,1,1,0,1,1,0,1,0,1,1,0,1,0,1,0,0,0,0,1,1,1,1,0,1,0,1,1,0,1,0,0,1,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,1,1,0,0,1,0,0,0,1,1,0,0,0,0,0,1,1,1,1,0,0,1,0,1,0,0,1,1,1,1,1,1,0,0,1,0,1,1,0,1,1,1,1,1,1,0,0,1,1,0,0,0,1,0,1,0,1,1,0,1,1,1,0,1,1,0,0,0,1,1,0,1,1,0,0,0,0,1,1,1,1],
[1,0,1,1,0,0,0,0,0,0,1,1,0,1,1,0,0,1,0,0,0,0,0,0,1,1,1,0,1,1,0,0,0,0,1,1,0,1,0,1,0,1,1,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1,1,0,0,0,1,0,1,1,1,1,0,0,0,0,1,0,1,0,0,0,0,1,1,1,1,1,1,1,0,0,1,0,0,0,0,1,1,0,1,1,1,1,0,0,1,0,1,1,0,0,1,0,1,0,1,0,1,1,1,1,0,0,1,1,0,0,0,1,1,0,0,0,1,0,1,0,1],
[1,1,0,0,0,1,0,1,1,0,0,1,1,0,0,0,1,1,0,1,1,1,0,0,1,0,0,1,1,1,1,0,0,0,1,0,0,1,0,0,1,0,0,0,0,1,1,1,0,0,1,1,0,0,1,1,1,0,1,0,1,1,0,0,0,1,1,0,1,1,0,0,0,1,0,1,0,1,1,1,0,0,0,1,0,1,1,0,0,1,1,0,1,1,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,1,0,1,1,1,1,0,1,1,0,1,0,0,1,0,0,0,1,1,1,0,0,1,1,1,1],
[0,0,0,0,1,0,0,1,1,0,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,0,1,1,0,0,0,1,1,0,0,1,0,0,0,1,0,0,1,1,1,1,0,0,1,1,0,0,0,1,0,1,1,0,0,0,0,0,0,1,0,0,1,1,1,0,0,1,0,1,1,0,0,0,1,1,1,1,0,0,1,1,1,1,0,0,1,1,0,0,0,1,1,0,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,1,1,1,0,1,0,1,1,1],
[1,0,0,1,1,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0,1,1,1,1,1,1,1,1,0,1,0,1,0,1,1,1,0,0,0,1,1,0,0,1,1,0,0,1,0,0,1,1,0,0,0,0,0,0,1,1,0,1,1,0,1,0,1,0,1,0,1,0,1,1,1,0,1,0,0,1,1,0,0,0,0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,1,0,0,0,1,1,1,1,1,1,0,0,0,1,1,1,1,1,0,0,0,1,1,1,0,0,0,1,0,0,1,0,1],
[0,1,1,0,0,1,1,0,1,1,0,1,1,1,1,1,1,1,1,0,0,1,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,1,0,0,1,0,0,0,0,1,1,1,0,0,1,0,0,1,0,0,1,1,1,0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,0,1,1,1,1,0,1,1,1,0,1,1,1,1,0,0,1,1,1,1,1,0,1,1,1,0,1,1,0,1,0,1,1,0,0,1,0,0,0,1,0,0,1,1,1,0,1,0,0,1,1,1,1,0,1,0,0,0,1,0,0],
[1,1,0,1,0,0,1,0,0,0,1,0,1,1,0,0,1,0,0,0,1,0,1,1,1,1,1,1,0,1,0,1,1,1,0,0,1,1,1,1,0,0,0,0,0,1,0,0,1,1,1,0,1,0,0,1,0,1,1,0,0,0,1,1,1,0,1,1,0,1,1,1,0,0,1,1,1,0,1,0,0,1,0,1,1,0,0,0,0,0,1,0,0,1,1,1,1,0,0,1,1,0,1,1,0,1,1,0,1,0,1,0,0,1,0,1,1,0,1,1,0,1,1,1,0,1,1,1,1,1,1,0,0,1,0],
[0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,1,1,0,1,1,0,1,0,1,0,0,0,1,1,1,0,0,0,1,0,0,1,0,0,0,1,0,0,0,1,1,1,1,0,1,0,1,1,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,1,0,1,1,0,0,0,1,0,0,0,0,0,1,0,1,1,0,0,1,1,1,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0],
[0,0,0,1,0,1,1,0,1,1,0,0,1,0,0,1,1,0,0,1,0,1,0,1,1,0,1,0,1,0,0,1,1,0,1,1,0,0,0,1,0,0,1,1,1,1,0,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,1,1,1,0,1,1,1,0,0,1,1,1,0,1,0,1,0,0,1,0,1,0,0,1,1,0,0,0,0,1,0,0,1,1,1,0,1,0,1,1,0,1,1,1,0,1,1,1,1,1,1,0,0,1,0,0,0,1,1,1,1,1,0,1,1,0,0,1,1,0,0],
[0,0,0,1,1,0,0,1,1,1,1,0,1,1,1,1,0,0,0,1,0,0,0,0,1,0,1,1,1,1,0,1,1,0,1,0,0,0,1,1,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,0,0,0,0,1,0,1,1,0,1,0,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,1,0,1,0,0,0,0,0,1,1,1,1,1,0,1,0,0,1,1,0,0,1,1,1,0,1,1,1,0,1,1,1,1,0,0,1,0,1,0,1,0,0,0,1,0,0,0],
[1,1,0,1,1,0,1,1,1,0,0,0,1,1,1,1,0,0,1,1,0,0,1,0,1,1,1,0,1,1,0,0,0,0,1,1,0,0,0,1,1,0,0,0,0,1,1,0,1,1,1,0,1,1,1,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,0,0,1,1,1,1,0,0,0,0,1,0,0,1,1,0,1,0,0,0,1,1,0,0,1,1,1,1],
[1,1,1,0,0,1,0,1,0,0,0,0,1,1,0,1,1,0,1,0,0,1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,0,1,1,0,0,0,1,1,1,0,1,1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,0,1,0,0,1,0,1,1,1,0,1,0,1,1,1,1,1,1,0,0,1,0,0,0,0,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,1,0,1,0,0,1,1,1,0,0],
[1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,1,1,0,1,1,0,0,0,0,1,0,0,0,1,0,0,1,1,1,1,0,1,1,1,1,1,0,1,0,1,1,0,0,1,0,0,0,0,1,1,0,0,1,1,0,0,1,1,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,1,0,0,1,1,1,1,0,1,0,1,0,0,1,0,1,1,1,1,1,0,1,0,0,0,1,0,1,1,0,1,1,1,1,0],
[1,1,0,1,1,1,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,1,0,1,0,1,0,1,0,0,1,0,0,1,1,1,1,1,1,0,0,0,0,1,1,0,1,0,1,1,0,1,1,1,1,0,1,0,0,1,1,1,1,1,0,0,0,0,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,0,0,1,1,0,1,1,0,0,1,1,0,1,0,1,1,0,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,0,1,0,0,0,1,1,0,1,0,0,1,1,0,0,1,0,1],
[0,0,1,0,0,1,1,1,0,0,1,1,0,1,1,0,1,1,1,1,1,0,1,0,1,1,0,0,0,0,0,1,1,0,1,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,0,0,1,1,0,0,0,0,1,1,1,0,1,0,0,0,0,1,0,1,1,1,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,1,0,0,0,1,1,0,0,1,1,0,1,0,0,0,1,0,1,1,1,1,1,0,1,1,1,0,1,0,0,1,0,1,1,1,1,1,0,0,0,1,0,1,0],
[0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,1,0,0,0,1,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,0,1,0,0,0,1,0,0,0,0,0,0,1,1,0,0,1,0,0,0,1,1,1,1,0,1,1,1,1,0,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,0,1,1,1,0,0,0,1,0,0,0,1,0,0,0,1,0,1,0,1,0,1,1,0,1,1,0,1,0,0,1,1,1,0,0,0,0,0,1,1,0,0,1,0,1,1,1,0,0,1,0],
[1,1,0,1,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,1,0,0,0,1,1,1,0,1,1,0,1,1,0,1,0,0,1,1,0,1,1,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,0,0,0,1,1,0,1,1,1,0,0,0,1,0,1,1,1,0,1,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,1,1,1,1,1,0,0,1,0,1,1,0,1,0],
[1,0,0,0,1,1,0,1,0,0,1,1,0,0,1,0,1,0,0,1,1,1,0,0,1,0,1,1,1,0,1,0,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,0,1,1,0,1,0,1,0,1,0,0,1,1,1,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0,0,1,1,1,1,1,0,1,0,0,1,0,1,1,1,0,1,0,1,1,1,1,0,0,0,0,1,0,1,0,1,1,0,1,0,1,0,1,0,1],
[1,1,0,0,0,0,1,1,1,1,1,0,1,0,0,1,0,0,1,0,1,0,1,1,0,1,0,1,1,1,1,1,1,1,1,0,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,0,0,1,1,0,0,1,1,1,0,1,1,0,1,1,1,1,1,1,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,0,0,1,1,1,0,0,1,0,1,0,1,1,1,1,1,0,0,1,1,0,1,1,1,1,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,0,0,1,1,1,0,1,0,0],
[0,0,1,0,1,0,1,0,0,1,1,1,1,1,1,0,0,1,1,0,1,0,1,1,1,1,0,0,0,1,1,1,0,0,0,1,0,1,1,1,1,1,0,1,0,1,0,1,1,1,1,0,0,1,0,0,0,1,1,0,0,1,1,1,0,1,1,0,1,1,1,1,1,0,1,0,0,0,1,1,0,0,1,1,0,0,0,0,1,1,0,0,1,0,1,1,0,1,1,1,1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,1,0,1,1,1,0,0,1,0,1,0,1,0,0,1,1,0,0,1,0],
[1,0,1,0,1,1,0,0,0,1,0,1,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,1,1,1,1,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,0,0,0,1,0,1,0,1,1,1,1,0,0,1,0,1,1,0,0,1,0,0,0,1,0,1,0,0,1,0,0,0,0,1,1,1,1,1,1,1,0,0,1,0,1,0,1,1,0,1,1,0,1,1,1,0,1,0,0,0,0,0,0,1,1,0,1,0,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0],
[0,1,1,0,1,1,1,0,1,0,1,0,0,0,1,0,1,1,0,1,0,1,1,1,0,1,1,0,1,0,1,0,1,1,1,1,1,1,1,0,0,1,0,0,1,1,1,1,1,1,0,0,0,0,1,0,1,1,0,1,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,0,1,0,0,0,0,1,1,0,1,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,1,1,1,1,0,0,0,1,0,0,0,0,1,1,1,1,0,1,0,0,0,1,1,0,1],
[0,1,1,0,1,1,0,1,0,1,0,1,0,1,0,0,1,1,1,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,1,0,0,1,1,0,0,1,0,0,1,0,0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,1,0,1,1,1,0,0,0,1,1,1,1,0,1,1,1,0,0,0,0,1,0,0,1,1,0,1,0,1,1,1,0,1,1,0,0,1,1,0,0,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,1,0,1,0,1,1,0,1,0,0,0,0,1],
[0,0,0,0,0,0,1,1,1,0,1,0,1,0,1,0,1,1,0,0,1,0,0,1,1,1,0,0,1,1,1,1,1,1,1,1,1,0,0,0,1,1,1,0,0,1,0,1,1,0,1,1,0,0,0,1,1,1,1,1,0,0,1,0,0,0,1,0,0,1,1,0,1,0,0,0,0,0,1,1,1,0,1,1,0,0,1,1,0,0,1,1,1,0,1,1,0,1,1,0,1,1,0,0,1,0,0,0,1,0,1,1,0,1,0,1,0,1,0,0,1,0,1,0,0,0,1,0,0,1,0,0,1,1,1],
[0,0,0,1,0,0,1,1,0,1,0,1,0,1,0,0,0,0,1,0,1,0,0,1,0,0,1,0,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,0,1,1,0,0,0,0,1,1,1,1,0,0,1,1,1,0,1,0,0,0,0,1,0,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,0,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,1,1,0,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,0,0,0,1,0,1,0,1,0,0,1,1,1,1,1,1,0],
[1,0,0,1,1,0,0,1,0,1,1,0,0,1,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,0,0,1,0,1,0,1,1,1,1,0,0,1,0,1,1,0,0,1,0,0,0,1,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,0,0,1,1,1,1,1,0,1,1,0,1,0,1,1,1,0,1,0,1,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,1,0,1,1,0,0,0,0,1,0,0,1,0,0,0,1,0,1,1,0,1,0,0,1,0],
[1,1,0,0,1,1,0,0,0,1,0,1,0,0,1,1,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,0,0,0,0,0,1,0,1,1,0,0,1,1,0,0,1,1,0,0,0,1,1,0,0,1,1,0,1,0,0,0,0,1,0,0,1,1,1,1,1,1,1,0,1,0,1,0,0,0,1,0,1,1,1,0,1,1,0,1,0,1,1,0,0,0,0,1,1,1,1,1,0,0,1,1,1,0,0,1,0,1,0,1,0,1,1,1,1,1,0,0,0,1,0,0,0,0,1,0,0,1,0],
[0,0,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,1,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,0,1,0,1,0,0,1,1,1,0,1,1,0,1,0,0,1,0,1,0,0,1,1,1,0,0,0,1,0,1,0,0,1,1,1,0,0,0,1,0,1,0,1,1,0,1,0,0,1,0,0,1,0,0,1,1,1,1,0,0,1,1,1,0,0,1,0,0,1,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,0,0,0],
[1,0,1,1,1,0,1,1,1,1,0,1,0,1,1,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,1,0,0,0,1,0,0,1,1,1,1,1,1,1,1,1,0,0,1,0,0,1,1,0,1,1,1,0,1,1,0,1,0,1,0,1,0,0,1,1,0,1,0,1,0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,1,0,1,1,0,1,0,1,0,1,1,1,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,1,0,1,1,1,1,1,0,0,0,1,1,1,0,1,0],
[0,0,1,0,1,1,1,1,0,1,1,0,1,1,1,1,0,0,1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,1,1,1,0,0,1,0,0,1,0,0,0,0,0,1,1,1,0,1,0,1,0,1,1,0,0,1,0,1,1,0,0,1,0,0,0,1,0,0,0,0,1,1,1,0,1,1,1,0,0,0,0,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,1,0,0,1,1,0,0,1,1,0,0,0,1,1],
[1,0,1,1,1,1,0,1,1,0,1,0,0,1,1,1,0,0,1,0,0,1,0,1,0,0,1,0,1,0,0,1,1,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,0,0,0,1,1,0,1,1,1,1,0,0,0,0,1,0,0,1,0,1,0,1,1,1,0,0,0,1,0,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,1,1,1,1,1,0,0,0,1,0,0,1,0,1,1,1,0,1,0,0,1,0,0,1,1,1,0,0,0,0,0,1,0,0,0,1,1,1,0,1,0,1,1],
[0,0,0,1,1,1,0,1,1,0,0,0,1,0,1,0,0,1,0,0,0,1,0,0,0,1,1,0,0,1,0,1,1,1,1,0,0,1,0,1,0,0,1,0,0,0,1,0,1,0,1,1,0,1,0,1,0,0,0,0,0,1,1,0,0,1,1,0,1,1,1,0,1,1,0,0,0,1,0,1,0,0,1,0,0,0,1,1,1,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,0,0,0,1,0,0,1,0,1,1,1,0,0,1,0,0,0,1,0,1,0,0,1,0,1,0,0,1,1,1],
[1,0,0,1,1,0,0,0,1,0,1,1,1,1,1,1,1,0,1,1,0,0,0,0,1,1,1,1,1,1,0,1,1,1,0,1,0,0,0,0,1,1,0,0,0,1,1,1,1,0,1,0,0,0,0,1,0,0,0,1,1,1,1,1,1,0,0,1,1,0,0,0,1,1,1,0,0,0,0,0,1,1,1,0,1,1,1,0,0,0,1,1,0,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0,0,1,0,1,1,0,1,1,0,0,1,0,0,1,1,1,1,1,1,1,0,0,0,1,0,0,0],
[1,1,1,0,1,0,1,0,1,1,0,1,1,1,1,0,1,1,0,1,1,1,0,0,0,0,1,1,0,0,1,1,1,1,0,1,0,0,1,0,1,1,1,0,1,0,1,0,1,0,0,0,0,0,0,1,1,0,0,1,1,1,0,0,1,0,0,0,1,0,0,1,1,1,1,1,1,0,0,1,0,0,1,1,1,1,1,0,1,0,1,0,1,0,0,0,0,1,1,1,1,1,1,1,0,0,1,0,0,1,1,0,1,1,0,0,1,0,0,1,1,1,0,1,1,0,1,1,0,0,1,0,1,1,1],
[1,1,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,1,0,0,1,1,1,1,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,1,0,0,0,0,1,1,1,1,1,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,0,0,1,0,0,1,1,1,0,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0,1,0,0,0,1,0,0,1,1,0,1,0,1,0,0,0,0,1,0,0,0],
[0,1,1,0,1,0,1,0,0,1,0,1,0,0,1,1,0,0,0,1,0,1,1,1,1,1,1,1,0,0,1,0,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,0,0,1,1,0,1,1,0,0,0,0,1,0,1,0,1,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,1,0,0,1,0,1,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,0,1,0,1,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0],
[1,1,1,0,1,1,0,1,1,0,0,1,0,1,1,1,0,0,1,1,0,1,0,0,0,1,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,1,1,1,0,1,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,1,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,0,1,0,1,1,1,0,1,1,0,1,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,0,0,1],
[0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,1,0,0,1,1,1,1,0,1,0,1,0,0,0,0,1,1,0,1,0,1,1,1,0,0,0,0,1,0,1,1,0,0,1,1,1,0,1,1,0,1,1,1,1,0,1,0,0,1,0,1,0,0,0,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0,1,0,0,1,1,1,0,1,1,0,0,0,1,1,0,0,0,1,1,1,1,1,0,0,0,0,0,1,0,1,1,1,0,0,0,1,0,1,1,1,0,0,0,1,0,0,0,1,0],
[1,1,0,0,0,1,0,0,0,1,0,1,1,1,0,0,1,1,0,1,1,0,1,0,0,0,0,0,1,1,0,1,1,1,1,0,1,1,1,0,0,0,1,1,0,0,1,0,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1,0,1,0,0,1,1,1,0,0,1,1,1,0,0,1,0,0,0,1,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1,1,0,0,1,1,0,0,0,1,0,1,1,1,0,1,1,0,1,0,1,0,0,0,0,1,0,0,1,0,0,1,0,1,1,1,1,1,0],
[0,0,1,1,1,0,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,0,1,1,1,0,0,1,0,1,0,0,1,0,0,0,1,1,1,0,1,0,0,1,1,1,0,1,1,0,0,1,0,0,1,0,1,1,0,0,0,1,1,0,1,1,0,1,0,1,1,1,0,1,0,1,1,0,0,0,1,0,1,1,1,1,1,1,0,1,1,1,0,0,0,0,1,0,1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,0,0,1,0,1,1,1,0,1,0,0,0,1,1,1,0,1,0,1,1,0,0],
[1,1,1,1,0,1,1,1,0,1,0,1,0,0,1,0,1,0,0,0,0,1,1,1,0,0,1,0,0,1,0,1,1,0,1,1,1,0,1,0,1,1,1,0,1,1,0,0,0,1,0,0,0,0,1,1,0,1,1,0,1,0,0,1,1,0,1,1,0,1,1,0,0,1,0,1,1,1,1,0,1,1,0,0,1,1,0,1,1,0,0,0,1,0,1,0,0,1,0,0,1,1,0,1,1,1,1,0,0,0,1,1,1,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,1,0,0,0,1,0],
[0,0,1,1,1,1,1,1,1,0,0,0,1,1,1,1,0,1,0,0,1,1,1,1,1,0,1,1,0,1,0,1,1,1,1,1,1,1,1,1,1,0,0,0,1,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,0,1,0,0,1,0,0,1,1,0,1,1,0,0,0,0,0,1,0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1,0,0,1,1,1,1,1,1,1,0,0,0,1,0,0,1,0,0,1,1,1,0,0,1,1],
[1,0,1,0,1,1,0,1,0,0,0,0,0,1,0,0,0,0,1,0,1,0,0,0,1,1,0,1,1,0,0,0,1,0,1,0,1,1,1,0,1,1,1,0,0,1,1,0,0,0,1,0,1,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,1,1,0,0,1,1,1,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,1,1,1,0,0,1,0,1,0,0,0,1,0,1,1,1,1,1,0,0,0,0,1,1,0,1,1,1,0,0,0,1,0],
[1,1,1,1,1,1,0,0,0,1,1,0,0,0,1,0,0,0,1,1,1,1,1,1,1,1,1,0,0,1,1,0,1,1,1,1,1,0,0,1,0,1,1,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,1,1,1,0,0,0,1,0,1,0,0,1,0,0,0,1,0,0,0,1,1,0,1,0,1,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,1,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,0,0,1,1,0,1,1,0,0,0,1],
[0,1,1,0,0,0,0,0,1,1,1,1,0,0,1,1,1,1,1,0,1,1,1,1,1,0,0,0,1,0,0,0,0,1,0,1,0,0,1,1,1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,1,1,1,0,1,0,1,0,1,0,1,1,0,1,1,0,1,1,0,0,0,0,1,0,0,1,1,0,0,1,1,0,1,1,1,0,1,0,1,0,1,0,1,1,1,1,1,1,1,0,0,0,1,0,1,0,1,1,1,1,0,1,1,0,1,1,0,0,1,1,1,1,0,1,0,1,1],
[0,1,1,0,0,1,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,1,0,0,0,1,1,1,0,1,1,0,0,0,0,1,0,0,1,0,0,1,1,1,0,0,0,1,1,1,1,0,1,0,1,0,1,1,0,0,0,1,1,0,0,1,1,1,0,0,0,0,0,1,1,0,1,0,1,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,1,0,0,1,0,1,0,1,0,0,1,1,1,0,0,0,1,1,1,1,1,1,0,1,1,0,0,1,0,1,1,0,0,1,0,1,0,1],
[1,0,1,1,1,1,0,0,1,0,1,0,0,1,1,1,0,0,1,1,0,1,0,0,1,0,1,1,1,0,1,0,0,1,0,1,1,0,1,0,1,1,1,0,0,0,0,1,1,1,1,1,1,1,0,1,1,1,1,0,0,1,1,0,1,0,1,1,1,0,1,1,0,0,0,0,0,1,1,1,0,0,0,0,1,0,1,1,1,0,0,1,1,1,1,0,1,1,1,1,0,1,1,0,1,0,1,0,0,0,1,0,1,0,1,1,1,1,0,1,0,1,0,0,1,1,0,0,1,1,0,0,1,1,1],
[1,1,1,0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,1,0,1,1,1,1,1,1,1,0,1,0,1,0,0,0,0,0,1,1,1,1,1,1,1,1,0,1,1,0,0,1,0,1,1,0,1,1,1,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,0,1,0,1,1,1,1,1,0,0,1,1,0,0,1,0,0,1,1,0,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,1,0,0,1,1,0,1,0,1,1,1,0,1,1,1],
[1,1,1,1,0,0,0,1,1,1,0,0,0,0,1,0,0,1,1,0,0,1,1,0,1,0,0,1,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,1,1,0,1,0,1,0,0,0,1,1,1,0,1,1,1,1,0,1,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,1,0,1,1,1,0,1,1,1,1,0,0,0,1,0,0,1,0,1,1,1,1,1,0,0,0,0,1,0,1,0,0,1,1,0,1,0,0,1,0,0,1,1,0],
[0,0,1,0,1,1,1,1,0,0,1,1,1,1,0,1,0,0,1,1,1,1,1,0,1,0,1,0,1,0,1,0,0,1,0,1,0,0,1,0,0,1,0,1,1,0,1,0,0,0,0,1,1,0,0,1,0,1,0,1,1,0,0,1,0,0,0,1,0,0,1,1,0,0,1,1,0,0,1,1,1,1,0,1,0,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,1,1,0,1,0,0,1,1,0,0,0,1,1,1,0,0,1,1,1,0,1,0,1,0,0,1,0,1,0,1,1,0,1,0,0],
[1,0,0,0,1,1,0,1,1,0,1,0,1,1,0,0,0,1,1,0,1,1,0,0,1,0,0,1,1,0,1,0,0,1,1,1,1,0,0,0,1,1,1,1,1,0,1,1,1,0,0,1,0,0,1,0,0,0,1,1,1,0,0,1,0,0,1,0,1,0,0,0,1,1,0,0,1,0,0,1,0,1,1,1,1,1,0,1,1,0,1,1,0,0,0,1,1,0,0,0,1,1,0,0,0,0,1,0,1,0,1,1,0,0,1,1,0,0,1,0,0,0,0,1,1,1,1,0,0,1,1,1,0,0,1],
[1,1,1,1,1,1,0,1,1,0,1,0,0,0,0,0,0,1,1,0,0,1,1,1,0,1,0,0,1,1,1,0,0,0,1,1,1,1,0,1,1,1,0,0,1,0,1,0,1,1,0,0,0,1,0,0,1,1,0,1,1,0,0,1,1,0,1,0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,0,0,0,0,1,1,1,1,0,0,1,0,1,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,1,0,1,0,1,0,1,0,0,1,1,1,1,0,0,1,1,0,0,1,1,1,0,1],
[1,0,1,1,1,1,0,0,0,1,0,1,0,1,1,0,1,1,0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,1,0,0,0,0,1,1,0,1,1,0,0,1,0,0,0,0,1,1,0,1,0,1,0,0,1,0,1,1,1,0,0,1,1,0,1,1,0,0,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,0,0,0,1,0,0,1,0,1,0,0,1,1,0,0,0,0,1,1,1,1,1,0,0,0,0,1,1,1,0,1,1,1,0,0,1,0,0,0],
[1,1,1,1,1,0,0,1,1,1,1,1,0,1,1,1,1,0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,0,1,0,1,0,1,0,0,1,0,0,0,0,1,1,1,0,1,0,0,1,1,0,1,0,1,0,0,1,1,0,1,0,1,1,0,0,1,1,1,1,0,1,0,1,0,0,1,0,1,1,0,0,1,1,0,1,0,1,0,1,0,0,0,1,1,0,1,1,0,0,0,0,0,1,0,0,1,1,1,1,0,1,0,1,0,0,0],
[1,1,0,0,0,0,1,1,1,0,1,0,0,0,0,1,0,1,1,0,0,1,0,1,1,0,1,1,0,1,0,0,0,1,0,1,1,0,1,1,0,1,0,0,0,1,1,1,1,0,0,1,0,1,1,1,0,0,0,1,0,1,1,0,0,0,0,1,0,0,1,1,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,1,0,1,1,0,0,1,0,1,1,1,1,0,1,1,1,0,0,1,0,0,1,1,0,1,0,0,0,1,0,0,1,0,1,1,1,1,0,1,0,0],
[0,1,0,1,0,1,0,0,0,0,0,0,0,1,1,1,0,1,0,0,1,1,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,1,0,1,1,1,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,1,1,0,1,1,1,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,0,1,1,0,0,0,0,0,0,1,0,0,1,0,0,1,0,1,1,0,1,0,0,1,1,1,1,1,1,1,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,1,1,1,0],
[0,0,0,0,0,0,0,1,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,1,0,0,0,0,1,1,0,0,0,0,1,0,1,0,1,0,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,1,1,0,0,0,1,0,0,1,0,1,0,1,0,0,1,0,1,1,0,1,1,0,0,1,0,0,1,1,1,1,0,0,1,0,1,1,0,1,1,0,1,1,0,0,1,0,0,1,0,0,1,0,0,1,1,1,0,1,0,1,0,0,0,1,0,1,1,0,0,0,0,1,1,1,0,1,0,1,0],
[1,0,0,1,0,1,1,1,0,0,0,1,0,1,0,0,0,0,1,0,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,1,0,1,0,0,0,0,1,0,1,1,1,1,1,1,0,0,1,0,1,1,1,1,0,0,1,0,1,1,1,0,0,0,0,1,0,1,1,0,0,1,1,0,0,1,0,0,0,1,1,0,0,1,1,1,1,1,1,1,0,1,1,1,0,1,0,0,1,1,1,0,0,0,1,1,0,0,0,0,0,1,1,0,1,0,0,1,0,1,1,1,1,1,1,0,0,1,1,0],
[1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,0,1,0,1,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,1,1,1,0,0,0,0,1,1,1,1,1,0,1,1,1,0,0,1,1,1,1,1,0,0,1,0,0,1,1,1,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,1,0,0,0,0,0,0,1,1,1,0,1,1,1,0,1,1,0,0,1,1,1,0,0,1,0,1,1,0,1,1,1],
[1,0,0,1,0,0,0,0,0,0,1,1,0,1,0,1,0,1,0,0,0,1,0,1,1,0,1,0,1,0,0,0,0,0,1,1,1,0,1,0,1,1,0,1,0,0,0,1,1,0,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,0,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,1,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,0,0,1,0,0,0,1,0,1,1,0,1,0,1,0,1,1],
[0,0,1,0,0,1,0,1,1,0,0,1,1,1,0,1,0,1,1,0,0,0,1,1,1,1,1,0,0,0,0,1,1,0,0,1,1,0,1,0,0,1,0,0,1,1,1,1,0,0,1,0,1,1,0,0,0,0,1,1,1,1,1,0,0,1,1,0,0,1,1,0,1,0,0,1,1,0,0,0,1,1,1,0,1,0,1,0,1,0,0,0,1,0,1,1,0,0,0,0,1,0,0,0,1,1,1,1,0,1,0,0,1,0,0,1,0,0,1,0,0,0,1,0,0,1,1,1,1,0,0,1,0,1,0],
[1,1,1,0,0,1,1,0,0,1,0,1,1,0,0,0,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0,1,0,1,0,1,1,1,0,0,0,1,1,0,1,1,1,1,0,0,1,0,1,0,0,0,0,1,0,1,0,1,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0,1,1,0,0,0,1,0,1,0,0,1,0,0,0,0,1,0,0,0,1,1,1,1,1,0,1,1,1,1,1,0,1,1,0,1,1,0,1,0,0,1,1,1,1,0,0,0,1,1,0,0,1,0],
[1,1,0,1,0,1,0,0,1,1,0,1,0,1,0,1,1,1,1,1,0,1,0,0,1,0,0,1,1,1,0,1,0,1,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,0,1,0,1,1,1,1,0,1,1,1,0,0,0,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,1,0,1,0,1,0,0,1,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,1,1,0,0,0,1,1],
[0,1,0,0,1,0,0,0,0,0,1,0,1,0,1,1,0,1,0,1,1,1,1,0,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,0,0,1,1,1,0,0,1,1,1,0,1,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,1,1,1,0,0,1,1,0,1,0,1,0,0,1,1,1,1,0,1,0,1,1,0,1,1,0,1,0,1,0,0,0,1,0,1,1,1,0,0,0,1,0,1,0,0,0,0,1,0,1,1,1,1,0,1,0,1,0,0,0,1,0,1],
[1,0,0,1,1,0,0,0,1,0,0,1,0,1,0,0,1,1,1,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,0,0,1,0,0,1,1,1,1,1,0,0,0,1,0,0,1,1,1,0,0,1,0,1,0,1,0,1,1,0,0,0,0,0,0,0,1,0,1,1,1,0,0,0,1,1,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,1,0,0,0,0,1,1,1,1,1,0,1,1,1,0],
[0,1,1,1,1,1,1,0,0,1,0,1,1,1,1,0,0,0,1,1,0,0,0,0,1,1,0,0,1,0,0,0,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,0,1,0,0,1,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,0,1,0,0,0,1,1,1,1,1,0,1,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,1,0,1,0,0,0,1,0,1,1,1,1,0,1,1,0,1,0,0,1],
[1,0,1,1,1,1,1,1,0,0,0,0,0,1,0,0,1,1,0,1,0,0,1,0,0,1,1,0,1,1,1,1,0,0,0,1,0,0,1,0,1,0,1,1,0,0,1,1,0,0,1,0,0,1,1,0,1,0,1,1,0,1,1,0,0,0,1,0,1,0,1,1,1,0,0,1,0,1,1,0,0,0,1,1,0,1,1,0,1,0,1,0,1,1,0,0,1,1,1,0,1,1,0,0,0,0,0,1,0,1,0,1,0,1,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,0,0,0,1,1,0],
[0,1,0,1,1,0,0,1,1,1,1,1,0,0,1,1,0,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0,0,1,1,0,0,1,0,1,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,0,1,0,0,1,0,1,0,1,1,0,0,0,1,0,1,0,1,1,1,0,0,0,0,0,1,0,1,0,1,1,1,1,0,1,1,1,1,0,0,0,1,1,0,1,1,0,1,1,0,1,0,1,1,0,0],
[1,1,1,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,0,1,1,1,1,0,0,1,1,0,1,1,1,1,0,1,0,0,1,0,1,1,1,1,0,1,1,0,1,1,1,1,1,0,1,0,1,1,0,1,1,0,0,1,0,1,0,1,0,1,1,1,0,1,1,0,0,0,1,1],
[0,1,0,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,1,1,0,1,1,1,1,0,0,1,1,1,1,1,1,0,1,0,1,1,0,0,1,0,1,1,0,0,1,1,1,1,0,1,0,0,0,1,0,1,0,0,1,1,0,0,0,1,0,0,0,0,1,0,1,1,0,1,0,0,0,1,1,0,0,0,1,0,1,0,1,1,1,0,1,1,0,0,0,1,0,1,1,0,0,1,0,0,1,1,1,0,0,1,0,0,0,1,0,0,1,0,1,0,0,0,1],
[0,1,0,1,1,0,1,1,1,1,1,0,0,1,0,1,0,1,0,1,1,1,0,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,1,0,0,1,1,0,0,0,0,1,0,0,1,0,1,0,1,0,1,0,0,1,0,1,0,0,0,1,1,1,0,0,0,1,1,1,0,1,0,1,1,0,1,0,1,1,0,0,1,0,0,1,1,1,1,0,0,0,1,0,1,0,1,0,0,1,0,0,1,0,1,1,1,1,1,0,1,1,0,0,0,0,0,1,1,0,1,1,0,0,1,1,0,0,1,0,0],
[0,0,0,1,1,1,0,0,1,0,0,0,1,1,0,1,1,1,1,1,1,1,0,1,0,1,1,1,0,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0,1,1,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,1,0,0,0,1,0,1,1,1,1,1,0,1,1,1,0,0,0,0,0,1,0,1,0,1,0,1,0,0,1,1,1,0,1,0,0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0,0,1,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0,1,0],
[0,1,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,1,1,0,1,1,1,1,0,0,1,0,0,1,0,1,0,1,0,1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,1,1,1,0,1,1,0,0,0,0,0,1,0,1,1,0,0,0,0,1,1,1,1,0,1,1,1,0,0,0,1,1,1,1,1,0,1,1,0,1,0,0,1,0,1,1,0,1,1,1,0,0,1,0,1,1,0,0,0,1,0,0,1,1,1,0,0,1,1,0,1,0,0],
[1,0,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,0,1,1,1,0,0,1,0,0,1,1,1,1,0,0,1,0,1,0,1,0,0,0,0,1,0,0,0,1,1,0,1,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,1,1,1,1,0,0,1,0,1,1,0,1,1,0,1,0,0,0,0,0,0,1,0,0,1,0,1,0,1,0,0,1,1,0,1,0,1,1,0,0,1,1,0,0,1,0,0,0,0,0,1,0,0,1,1,0,0,1,0,1,1,0,1],
[0,1,0,1,1,0,0,0,0,1,1,0,0,1,1,0,1,0,1,1,0,1,1,0,0,1,1,1,1,0,1,1,1,1,1,1,1,0,1,0,1,1,1,0,0,1,1,1,1,0,0,0,0,1,1,1,0,0,0,1,1,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0,0,1,1,0,0,1,1,0,0,1,0,0,0,1,1,1,1,0,1,1,1,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0],
[0,1,0,0,0,0,0,1,1,1,0,1,0,1,1,0,0,0,0,1,1,0,1,1,0,0,0,0,1,1,1,1,0,0,1,1,0,1,0,0,0,1,1,0,1,0,0,0,0,1,1,0,0,1,1,0,0,1,0,0,0,1,1,0,0,1,1,1,0,1,0,0,0,0,0,0,1,1,1,0,0,1,1,0,1,1,1,0,0,0,0,1,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,1,0,1,1,0,0,1,0,0,0,0,0,1,1,1,0,1,0],
[0,1,0,1,0,0,1,1,0,1,1,0,1,1,0,1,0,0,1,0,0,1,1,0,0,1,1,1,0,1,0,1,1,1,0,0,1,1,0,1,0,0,1,1,1,1,0,0,0,0,0,0,1,1,0,0,1,0,1,0,0,1,0,0,0,1,0,0,1,0,1,1,0,1,0,1,0,0,0,0,0,0,0,0,1,1,1,0,0,1,1,1,1,0,0,0,1,0,0,0,0,1,0,0,1,1,1,0,0,1,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,0,1,1,1,1,0,0,0,1],
[0,0,1,0,1,1,0,1,1,0,0,0,1,0,0,0,0,1,0,0,1,1,1,1,1,1,1,0,1,1,0,1,0,1,1,0,0,0,0,0,0,1,0,1,1,1,1,1,0,0,0,0,0,1,0,1,1,1,1,0,0,1,1,0,0,0,0,1,1,1,0,0,0,1,1,1,0,0,0,1,0,0,0,1,1,1,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,1,1,0,1,1,0,1,1,1,0,0,1,1,1,1,0,1,0,0,1,0,0,0,1,0,0,0,0,1,1],
[1,1,1,1,1,1,0,0,0,0,1,1,1,0,0,0,0,1,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,0,1,0,1,1,1,0,1,0,0,0,1,1,1,1,1,1,0,1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,1,1,0,0,1,0,0,1,0,1,1,1,1,1,0,1,1,0,0,1,0,0,0,0,1,0,0,0,1,0,1,0,0,1,0,1,1,1,0,0,0,1,1,1,0,0,0,0,1,0,1,1,1,0,1,0,0,0,1,1,1,1,0,1,1,1,1,1],
[1,1,0,1,1,1,1,1,0,1,0,0,1,0,0,1,0,0,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,1,0,0,0,0,0,1,1,1,1,1,0,0,0,1,1,0,1,0,1,1,1,0,0,1,1,0,1,0,1,0,0,1,1,1,1,1,0,1,1,0,1,0,0,1,1,1,1,0,1,0,0,0,0,0,1,1,0,0,0,1,0,1,1,0,1,1,1,0,1,1,1,1,1,0,1,1,0,1,0,0,0,0,0,0,1,1,1,0,1,0,1,1,0,0,0,1,1,0,0,0],
[0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,0,0,1,1,0,1,1,1,1,0,1,1,0,0,0,1,0,1,0,1,1,0,0,1,1,0,1,1,0,0,0,1,0,1,1,1,1,0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,1,0,0,1,1,1,0,1,0,1,1,1,1,1,0,0,1,0,1,1,0,0,1,1,1,1,0,0,0,0,0,0,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,1,0,1,1,0,0],
[0,0,0,1,1,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,1,1,0,0,1,0,0,0,0,0,0,1,0,1,0,1,0,0,1,0,1,1,0,1,0,1,1,0,0,1,0,1,0,0,1,0,1,0,1,1,0,1,0,0,1,1,0,0,1,1,0,0,0,0,1,0,0,1,0,0,0,0,1,0,1,0,1,0,0,0,1,1,0,1,0,0,1,1,0,0,1,0,0,1,0,0,1,1,1,1,0,1,1,0,0,0,1,0,0,1,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0],
[1,0,0,0,0,1,0,1,1,1,1,0,1,0,1,0,1,1,0,1,0,0,0,1,1,1,0,0,0,0,1,0,0,1,0,1,1,1,0,0,1,1,0,0,0,0,0,1,1,1,1,0,1,0,0,0,0,1,0,0,0,0,0,1,0,1,1,1,0,0,1,1,1,0,0,0,0,1,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,0,1,1,1,1,0,0,0,1,0,0,1,1,0,1,0,0,1,0,1,1,0,0,0,1,0,1,1,0,0,0,1],
[0,1,0,1,1,1,0,0,1,0,0,0,1,0,1,0,1,0,0,1,1,1,1,1,0,0,0,1,1,1,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,0,0,0,1,1,0,1,0,1,0,1,1,0,0,1,0,0,1,0,1,0,1,0,1,1,1,1,1,0,0,0,1,0,1,0,0,1,1,0,1,0,0,1,1,1,1,1,1,1,0,1,0,0,0,1,0,1,1,0,0,0,1,0,0,1,0,1,1,0,0,0,1,0,0],
[0,1,1,0,0,0,1,1,1,1,0,1,0,0,0,0,0,1,1,0,0,1,0,1,1,1,1,0,1,1,1,1,1,1,1,0,1,1,1,0,1,0,0,1,0,1,0,0,1,0,1,1,0,0,1,1,1,0,1,0,1,1,1,0,0,1,1,0,0,1,0,1,1,0,0,0,0,0,0,1,1,1,0,0,1,1,0,1,0,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,1,1,1,0,1,0,1,0,1,1,1,1,1,0],
[1,0,1,0,0,0,0,1,1,1,1,1,0,0,1,1,1,1,0,0,1,0,0,0,1,0,1,1,0,0,1,0,1,0,1,1,1,0,1,0,0,0,0,0,1,1,1,0,1,0,0,1,1,0,1,0,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,1,1,1,0,0,1,0,1,1,0,0,1,1,1,0,1,0,1,0,1,0,0,1,0,1,1,0,1,0,1,1,0,0,0,0,1,0,0,1,1,1,1,1,0,0,1,0,1,1,0,1,1,1,0,0],
[1,0,1,1,0,0,0,1,1,1,1,0,1,0,0,0,0,0,1,0,1,1,1,0,0,0,1,0,0,0,1,1,0,0,0,1,1,1,1,0,1,1,1,0,1,0,0,1,1,0,0,1,1,0,0,1,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1,1,0,0,1,1,1,0,1,0,1,1,0,1,1,0,0,1,0,1,1,0,0,0,0,0,1,1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,1,1],
[0,0,1,0,1,1,1,0,0,0,1,0,0,1,1,1,1,1,0,0,0,0,0,1,0,0,0,1,0,1,1,1,0,1,1,0,0,0,1,0,1,0,0,0,1,1,0,1,0,0,0,0,1,0,0,1,0,1,0,1,1,0,1,1,1,0,1,0,1,1,0,1,1,1,0,0,1,0,1,0,1,0,0,1,0,0,0,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0,0,0,0,1,0,1,1,1,1,0,0,1,1,0,0,0,1,0,1,0,1,1,1,0,1,1,0,1,0,1,1,1,0],
[0,1,0,1,1,1,0,1,1,0,0,1,0,1,1,0,0,1,0,0,1,0,0,0,1,1,1,0,1,0,1,1,1,0,1,0,0,1,0,0,1,1,0,0,1,1,1,0,0,1,0,1,1,0,1,0,1,0,0,1,1,0,1,1,0,0,1,1,1,0,1,0,1,0,1,0,0,1,1,0,0,1,1,0,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0,1,1,1,1,1,1,0,0,1,1,1,1,0,1,1,0,0,1,1,0,1,1,0,1,1,0,0,1,0,1,0,1,0,0,1,0],
[1,1,0,1,0,0,1,0,0,0,0,1,1,1,0,1,0,0,1,1,0,0,1,0,1,1,1,0,1,1,0,0,1,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,1,0,0,0,0,1,1,0,0,0,1,1,0,1,0,1,1,1,0,1,0,0,1,1,0,0,0,1,0,0,0,1,1,0,0,1,0,0,0,0,1,1,1,0,1,1,1,0,1,0,0,0,0,1,1],
[0,0,0,0,1,0,0,0,1,0,0,1,0,0,0,1,1,0,1,1,0,0,1,0,0,0,0,1,0,1,0,1,0,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,0,1,0,1,1,0,1,1,0,1,1,0,1,1,0,0,1,0,1,0,1,0,1,1,0,1,0,0,0,1,0,1,1,1,0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,0,1,0,0,0,0,1,0,1,1,0,1,1,0,1,1,1,1,1,0,1,0,1,0,1,0,1,1,0,0,1,0,0,0,1,1],
[1,1,0,1,1,0,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,1,0,0,1,1,0,1,1,0,0,1,1,1,0,0,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,1,1,0,1,1,0,1,1,1,1,1,1,0,1,1,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,1,1,1,0,0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,1,1,1,0,0,1,1,0,1,1,0,1,0,1,1,1,0,1,0,1,0,1],
[1,1,1,1,0,0,1,1,1,1,0,1,1,1,1,0,1,1,0,0,0,1,1,1,0,1,1,0,1,1,1,0,0,1,1,0,0,1,1,1,0,0,1,0,0,1,0,0,0,1,0,1,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,0,1,0,0,0,0,0,0,1,0,0,0,1,0,1,1,1,0,1,0,0,1,1,1,0,1,1,0,1,0,1,1,1,1,1,1,0,0,0,0,1,0,0,1,1,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,0,1,1,0,1],
[0,1,0,0,1,0,0,0,1,1,1,0,1,1,0,0,1,0,0,0,1,1,0,0,1,1,1,1,0,1,0,0,0,1,0,0,0,1,1,1,1,0,0,0,1,1,1,0,1,1,1,1,0,0,1,1,0,0,1,1,0,1,1,1,0,1,0,0,1,1,1,0,0,1,1,1,1,0,1,1,0,0,0,1,1,0,0,1,0,1,1,1,1,1,0,1,0,0,0,1,0,1,1,0,0,1,0,0,0,1,0,1,1,1,0,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
[0,0,1,1,0,0,1,0,0,0,1,0,1,0,0,1,1,0,1,0,0,0,1,1,0,1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,1,0,1,1,1,0,0,0,1,1,1,1,1,1,0,1,1,1,1,0,1,0,0,0,1,1,1,1,1,1,0,1,0,0,0,1,0,0,1,1,1,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,0,1,0,0,0,0,1,1,0,0,1,1,0,1,1,0,1,1,0,0,1,1,1,1,0,1,0,0,0],
[1,0,0,1,1,0,1,1,1,1,1,1,0,1,0,1,0,0,1,1,0,1,1,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,0,0,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,0,0,1,1,0,1,0,0,0,1,1,1,0,0,1,1,0,0,0,0,1,1,1,0,0,1,0,0,1,1,1,1,1,1,1,1,0,0,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,0,1,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,1,0],
[0,0,1,1,1,0,1,0,1,1,0,1,0,0,1,0,0,0,0,1,0,0,1,0,1,0,1,1,0,0,0,0,0,0,1,1,0,1,1,1,1,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,1,0,1,1,0,0,1,0,1,1,0,1,1,1,1,1,1,1,1,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,1,0,0,0,1,1,0,0,0,0,1,1,1,1,0,1,1,0,0,1,0,0,0,0,1,1,1,0,1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0],
[0,1,0,1,1,1,0,1,1,1,0,1,1,1,1,0,0,0,0,1,1,0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,1,1,1,1,1,0,1,0,0,1,0,1,1,0,0,0,0,0,1,0,0,0,0,1,1,1,1,0,1,0,1,0,1,1,1,0,0,0,1,0,0,0,0,1,0,0,1,0,1,1,0,1,1,0,0,1,1,1,1,0,1,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,1,1,0,0,1,1,0,0,1,1,1,0,0,1],
[0,0,0,1,0,1,0,1,1,0,0,0,1,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,0,0,0,0,0,0,1,1,0,0,0,1,0,1,1,0,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,1,1,1,1,1,1,0,0,1,0,0,0,0,0,1,0,0,1,1,0,1,0,1,1,1,1,1,1,1,1,0,1,0,0,1,1,0,1,0],
[0,1,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,1,1,1,0,1,0,0,1,0,0,1,1,0,1,1,1,1,1,1,0,0,1,0,0,1,1,1,0,1,0,0,1,1,1,1,0,0,1,0,0,1,0,1,1,1,0,1,0,0,1,0,1,1,0,0,1,0,1,0,0,1,0,1,1,0,0,0,0,0,0,0,1,1,0,1,1,0,1,1,0,1,0,1,0,1,0,1,0,0,0,0,1,0,1,1,1,1,0,1,0,0,0,0,1,0,0,1,1,1,1,0,0,0,0,0,0,0,0],
[1,1,0,1,1,1,1,0,1,0,0,0,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,1,0,0,1,1,0,0,1,1,1,0,1,1,0,0,0,0,1,0,1,1,0,1,1,1,1,1,0,1,0,1,0,0,1,1,1,1,0,0,0,1,0,0,1,1,1,0,1,0,1,1,0,0,0,1,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,1,1,1,1,1,0,0,0,1,1,0,1,1,1,1,1,1,1,0,0],
[1,0,1,1,1,1,0,0,0,1,0,0,1,1,1,0,0,0,1,0,1,0,1,1,0,1,1,0,1,0,0,0,0,1,1,1,0,0,0,0,1,1,1,1,0,1,0,1,1,0,0,1,1,1,0,0,1,0,1,0,0,0,0,0,1,0,1,0,1,0,0,1,1,1,0,1,0,1,1,0,1,0,0,1,0,0,1,1,0,0,1,1,0,0,0,1,1,1,0,1,1,1,0,0,0,1,1,0,1,1,1,1,0,1,0,1,0,1,1,1,1,1,0,0,0,1,0,0,1,0,0,1,0,1,0],
[1,1,1,0,0,0,0,1,1,1,0,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,1,1,0,0,1,1,0,1,1,1,1,1,0,0,0,1,0,0,0,0,1,1,1,1,1,0,1,0,0,0,0,1,0,1,1,0,0,1,0,1,0,0,1,1,0,1,0,0,0,0,0,0,1,1,1,1,1,0,1,1,1,0,1,1,1,1,0,1,0,1,0,1,1,0,0,1,1,1,1,1,1,1,0,1,0,0,1,0,0,0,1,1,1,0,0,1,0,1,1,1,0,0,0,0,1,1,1,0,1],
[0,0,1,0,0,1,0,0,0,1,1,1,1,0,0,1,1,0,0,1,1,1,1,0,0,0,1,0,0,0,0,1,0,1,0,1,0,0,0,1,0,0,1,1,0,1,1,1,0,0,0,1,1,1,0,0,0,0,0,0,1,0,1,0,1,0,0,1,1,0,1,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,1,0,1,1,0,1,1,1,0,1,1,1,1,0,1,0,0,1,0,1,1,0,0,0,1,1,0,1,1,1,1,1,1,1,0,0,0,1,0,1,1,0,0,0,0,0,1,0]
]
# s^transpose
SD_270_0_s = vector(GF(2),[0,1,1,1,0,0,1,1,0,0,1,1,1,0,0,1,1,0,0,1,1,0,1,0,1,0,1,0,0,0,0,0,0,1,0,1,0,1,1,1,1,0,0,1,1,0,0,1,0,1,0,0,1,1,1,0,1,1,1,1,0,0,1,1,1,1,1,1,0,0,0,1,0,0,0,0,1,0,1,1,1,1,1,0,1,1,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,1,0,0,1,1,1,0,1,0,1,1,1,1,1,1,0,1,1,0,0,1,0,0,1])
Hm = matrix(GF(2),H)
Hm = Hm.transpose()
I = identity_matrix(GF(2),135)
H = block_matrix([I,Hm],nrows=1)
SD_270_0 = LinearCode(H)

