@RemoteFile(
    SBOMB_WATER, 
    "https://www.dropbox.com/scl/fi/7lvc5j0dqcbztcwd3rnnp/water.nc?rlkey=yhpmf5evtzrppr14v691l97j3&dl=1",
    file = "water.nc", 
    updates=:never,
    dir = joinpath(@__DIR__))
!isfile(SBOMB_WATER) && download(SBOMB_WATER)

@RemoteFile(
    SBOMB_WIND, 
    "https://www.dropbox.com/scl/fi/k3vrbuxegp6mj6y5oqpbo/wind.nc?rlkey=jgi2zvrwagc2e1tcwt1nvzz54&dl=0",
    file = "wind.nc", 
    updates=:never,
    dir = joinpath(@__DIR__))
!isfile(SBOMB_WIND) && download(SBOMB_WIND)

@RemoteFile(
    SBOMB_WAVES, 
    "https://www.dropbox.com/scl/fi/j5b6akflmgvhcfy124l4l/waves.nc?rlkey=me7c7tne70kd5xrt3oy7o15tm&dl=0",
    file = "waves.nc", 
    updates=:never,
    dir = joinpath(@__DIR__))
!isfile(SBOMB_WAVES) && download(SBOMB_WAVES)

@RemoteFile(
    SBOMB_TEMP, 
    "https://www.dropbox.com/scl/fi/oyrk1klf8rugk3mty0erx/currents-temperature.nc?rlkey=3t6znxletyret5jo0v8447dpf&dl=1",
    file = "currents-temperature.nc", 
    updates=:never,
    dir = joinpath(@__DIR__))
!isfile(SBOMB_TEMP) && download(SBOMB_TEMP)

@RemoteFile(
    SBOMB_NUTRIENTS, 
    "https://www.dropbox.com/scl/fi/pj3ovx8owlxqxaj43x2py/nutrients.nc?rlkey=jsulyltckanjqdmi8ygxrddwo&dl=1",
    file = "nutrients.nc", 
    updates=:never,
    dir = joinpath(@__DIR__))
!isfile(SBOMB_NUTRIENTS) && download(SBOMB_NUTRIENTS)