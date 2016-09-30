MDplot <- function(path1){
	fs=.Platform$file.sep
	count = 1
	setwd(path1)
	for(count in 1:3){
		if(count==1) {
			PathBonds <- paste(path1,"BONDS_plot",sep = fs)
			message(PathBonds)
			setwd(PathBonds)
		} else if(count==2) {
			PathAngles <- paste(path1,"ANGLES_plot",sep = fs)
			message(PathAngles)
			setwd(PathAngles)
		} else if(count==3) {
			PathD <- paste(path1,"DIHEDRAL_plot",sep = fs)
			message(PathD)
			setwd(PathD)
		} else {
			stop("Counter past 3")
		}
		for(i in list.files(path=".")){
			strin <- i
			message(i)
			title <- unlist(strsplit(strin, split="\\."))
			if(count==1) {
				fulltitle <- paste(title[1]," Bond Length")
				lab_y <-"Bond Length (A)"
				data <- try(read.table(i,sep=",",header=T),FALSE)
				attach(data)
				png(paste(title[1],".png"), width = 10, height = 7, units = 'in', res = 300)
				plot(Geom, Value, main=fulltitle, xlab="Time Step (1000's)", ylab=lab_y) 
				dev.off()
				detach(data)
			} else if(count==2) {
				fulltitle <- paste(title[1]," Bond Angle")
				lab_y <-"Bond Angle (degrees)"
				data <- try(read.table(i,sep=",",header=T),FALSE)
				attach(data)
				png(paste(title[1],".png"), width = 10, height = 7, units = 'in', res = 300)
				plot(Geom, Value, main=fulltitle, xlab="Time Step (1000's)", ylab=lab_y) 
				dev.off()
				detach(data)
			} else if(count==3) {
			fulltitle <- paste(title[1]," Dihedral Angle")
				lab_y <-"Dihedral Angle (degrees)"
				data <- try(read.table(i,sep=",",header=T),FALSE)
				attach(data)
				png(paste(title[1],".png"), width = 10, height = 7, units = 'in', res = 300)
				plot(Geom, Positive_sense_angle, main=fulltitle, xlab="Time Step (1000's)", ylab=lab_y) 
				dev.off()
				detach(data)
			}
			
		}
		message("Finished set")
		setwd(path1)
	}
	stop("FINISHED")
}
