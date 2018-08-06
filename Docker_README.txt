# Commands to build, run, start, stop and remove the container for PyBox

# To build a new image, run the following command in the directory of the supplied Dockerfile:
docker build -t pybox .
# Please note that when the build reaches the section on whether to install MS Visual Studio or not:
# 'Do you wish to proceed with the installation of Microsoft VSCode? [yes|no]', it can appear to hang. However this is seemingly a sympton # of the build process. Please do wait for completion.

# After this has completed [which may take some time], type the following to see your new image listed
docker images

# To create and run a new container based on this image, with a name 'project_pybox', type:
docker run --name=project_pybox -it pybox
# This will take you in to the container. So, lets run the gas phase model in PyBox whilst you are there.
# Change directory to where PyBox is located:
cd /Code/Git_repos/PyBox/
# Lets run the default simulation:
python Gas_simulation.py 
# Dont worry about the error message regarding the Matplotlib plots. This is a result of working in a Docker container
# lets leave the container. Type:
exit
# You are now back in your original 'host' operating system.

# To list all available containers, type:
docker ps -a
# You should see how long ago you exited our 'project_pybox' container. Lets restart this container.
# To do this type the following:
docker start project_pybox
docker exec -it project_pybox bash
# We are back in! Lets exit again:
exit

# To remove the container, you can type:
docker rm project_pybox
# To remove the image on which the container is based, which was built from our Dockerfile, type:
docker rmi pybox
 
