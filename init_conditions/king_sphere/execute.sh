#/usr/bin/bash

touch table_of_parameters.out
rm table_of_parameters.out # Supprime la table des données précédément créée

gcc king.c profil.c particules.c -o main.x -lm  # Compilation 

if [ $? != 1 ] # Si la compilation fonctionne
then

	echo " "
	echo -e "Compilation done"
	echo " "

echo "Exploring the input parameters or generation of particles ?"
select j in table particles; do # Faire un choix sur la table des paramètres ou le sortie avec N particules

        if [ x$j = xparticles ]; then # Si sortie avec N particules
				echo "Sigma = "	;read sigma
				echo "R0 = "	
				read ro
				echo "Ratio = "	
				read ratio
				echo "Nb_part ="
				read nb_part
				echo "Mu fonction (newton, simple, std, exp, pot_ext (g_ext (EFE) dominated only for simple) )?"
				read nu_func
				echo -e "$sigma">.instructions
				echo -e "$ro">>.instructions
				echo -e "$ratio">>.instructions
				echo -e "$nb_part">>.instructions
				echo -e "$nu_func">>.instructions
				if [ $nu_func = "pot_ext" ]; then
								echo "g_ext (km2.s-2.kpc-1) = "
								read g_ext
								echo -e "$g_ext">>.instructions
				fi				


				./main.x<.instructions>.out
				mv king.out output/
				mv ic_part output/
				echo " "
				more table_of_parameters.out
				echo " "
                break

        elif [ x$j = xtable ]; then # Si sortie = table des paramètres
			 echo -e "sigma\tro\tratio\t->\tC\tRt\tMtot">table_of_parameters.out
			for sigma in `seq 1 10`  # Boucle sur les trois paramètres d'entrée
			do
				for ro in `seq 1 20`
				do
				ro=`echo -e "($ro*0.1)" | bc`
					for ratio in `seq 1 10`
					do
						#ratio=`echo -e "($ratio*0.1)" | bc`
						echo -e "\rsigma = $sigma \tro= $ro \tratio= $ratio\c"
						echo -e "$sigma">.instructions
						echo -e "$ro">>.instructions
						echo -e "$ratio">>.instructions
						echo -e "10">>.instructions
						echo -e "simple">>.instructions
						echo -e "512.51">>.instructions
						./main.x<.instructions>.out
					done 
				done
			done
				mv table_of_parameters.out	output/			
				echo -e "\n"
                break

        else
			echo "Bad Awnser"
        fi
done


else # Si problème de compilation
	echo " "
	echo -e "Compilation failed"
	echo " "
fi
