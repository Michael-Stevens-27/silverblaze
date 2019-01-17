Bootstrap: debootstrap
OSversion: xenial
MirrorURL: http://archive.ubuntu.com/ubuntu/


%post
        mkdir /data
        apt-get update
        ## Add more repos
        apt-get install -y software-properties-common
        ## Add universe
        add-apt-repository universe

        apt-get update
        apt-get install -y vim git wget tar unzip sharutils sudo curl
        apt-get install -y libxml2-dev freeglut3-dev libgeos-dev libcurl4-openssl-dev libssl-dev libssh2-1-dev libfftw3-dev libfftw3-doc apt-transport-https libpng-dev
        add-apt-repository 'deb [arch=amd64,i386] https://cran.rstudio.com/bin/linux/ubuntu xenial/'
        apt-get update
        apt-get install -y --allow-unauthenticated r-base
	
	# download repo for GDAL
	add-apt-repository -y ppa:ubuntugis/ppa
	apt-get update
	apt-get install -y gdal-bin libproj-dev libgdal-dev

	# set locale
	locale-gen en_GB.UTF-8

        # R PACKAGES
        R --vanilla <<END
	Sys.setlocale("LC_ALL", "en_GB.UTF-8")
	        
	# CRAN packages
        Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

        install.packages("geosphere", repos='http://cran.us.r-project.org')
        install.packages("scales", repos='http://cran.us.r-project.org')
        install.packages("jpeg", repos='http://cran.us.r-project.org')
        install.packages("proto", repos='http://cran.us.r-project.org')
        install.packages("rjson", repos='http://cran.us.r-project.org')
        install.packages("devtools", repos='http://cran.us.r-project.org')
	install.packages("rgdal", repos='http://cran.us.r-project.org')
	install.packages("RgoogleMaps", repos='http://cran.us.r-project.org')

        library(devtools)  
	install_github("dkahle/ggmap")
        install_github("bobverity/Rgeoprofile")
	install_github("michael-stevens-27/silverblaze")
