Stage0 += baseimage(image='ubuntu:18.04')
  
Stage0 += packages(ospackages=['fortune', 'cowsay', 'lolcat'])
Stage0 += environment(variables={'PATH': '/usr/games:$PATH'})

Stage0 += label(metadata={'Author': '"Pawsey Supercomputing Centre"'})
Stage0 += label(metadata={'Version': 'v0.0.1'})

