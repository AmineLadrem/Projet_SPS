import sys
import time
sys.path.insert(0, './lib')

from rich.console import Console
from rich.markdown import Markdown
from rich.panel import Panel
from adn import *
import random

choice=None

def create_menu(options):
    console = Console()

    markdown_text = "# Centre de Commande\n"
    for i, option in enumerate(options, 1):
        markdown_text += f"{i}. {option}\n"

    markdown = Markdown(markdown_text)
    panel = Panel(markdown, border_style="blue", expand=False)

    console.print(panel)

    choice = int(input("Veuillez entrer le numero de votre choix: ")) - 1
    if 0 <= choice < len(options):
        console.clear()
        message = Markdown(f"## Vous avez choisi l'option : {options[choice]}")
        message_panel = Panel(message, border_style="green", expand=False)
        console.print(message_panel)
        return choice
    else:
        print("Invalid choice. Please enter a number between 1 and", len(options))
        return None
 
options = [
    "Générer une chaîne ADN aléatoire d’une longueur donnée",
    "Choisir une chaîne ADN à partir d’un fichier",
    "Afficher la séquence ADN",  # Added option to display ADN sequence
    "Vérifier la validité de la chaîne ADN (si lue à partir d’un fichier)",
    "Calculer les fréquences des bases nucléiques dans la chaîne ADN",
    "Transcrire la chaîne ADN en une chaîne ARN",
    "Transcrire la chaîne ARN résultante en protéines (acides aminés)",
    "Calculer le complément inverse de la chaîne ADN",
    "Calculer le taux de GC de la séquence ADN",
    "Calculer les fréquences de codons dans la chaîne ADN",
    "Réaliser des mutations aléatoires sur la chaîne ADN",
    "Chercher un motif dans la chaîne ADN",
    "Générer la chaîne ADN consensus et la matrice profil",
    "Quitter"
]

adn_sequence=None

def main():
    while True:
        choice = None
        while choice is None:
            choice = create_menu(options)

        if choice == 0:
            length = int(input("Veuillez entrer la longueur de la chaine ADN: "))
            adn_sequence = generate_adn(length)
            print(adn_sequence)
            input("Press Enter to go back to the menu...")
        
        if choice == 1:
            file_path = input("Veuillez entrer le chemin / nom du fichier: ")
            description, adn_sequence = read_adn_file(file_path)    
            print(description)
            print(adn_sequence)
            input("Press Enter to go back to the menu...")
            
        if choice == 2:
            print(adn_sequence)
            input("Press Enter to go back to the menu...")
            
        if choice == 3:
            is_adn_valid(adn_sequence)    
            input("Press Enter to go back to the menu...")
            
        if choice == 4:
            frequencies = bases_freq(adn_sequence)
            print(frequencies)
            input("Press Enter to go back to the menu...")
            
        if choice == 5:
            rna_sequence = translate_to_arn(adn_sequence)
            print(rna_sequence)
            input("Press Enter to go back to the menu...")
            
        if choice == 6:
            rna_sequence = translate_to_arn(adn_sequence)
            proteins = translate_to_proteins(rna_sequence)
            print(proteins)
            input("Press Enter to go back to the menu...")
            
        if choice == 7:
            reverse = reverse_complement(adn_sequence)
            print(reverse)
            input("Press Enter to go back to the menu...")
            
        if choice == 8:
            gc_content = calculate_gc_content(adn_sequence)
            print(gc_content)
            input("Press Enter to go back to the menu...")
            
        if choice == 9:
            codon_frequencies = calculate_codon_frequencies(adn_sequence)
            print(codon_frequencies)
            input("Press Enter to go back to the menu...")
            

        if choice == 10:
            mutation_rate = float(input("Veuillez entrer le taux de mutation: "))
            mutated_sequence = mutate_adn_sequence(adn_sequence, mutation_rate)
            print(mutated_sequence)
            input("Press Enter to go back to the menu...")
            
            
        if choice == 11:
            motif = input("Veuillez entrer le motif à rechercher: ")
            positions = adn_motif(adn_sequence, motif)
            print(positions)
            input("Press Enter to go back to the menu...")
            
        if choice == 12:
            consensus_sequence, profile_matrix = generate_consensus_profile(adn_sequence)
            print(consensus_sequence)
            print(profile_matrix)
            input("Press Enter to go back to the menu...")
            
        if choice == 13:
            print(adn_sequence)
            input("Press Enter to go back to the menu...")
            
        if choice == 14:
            input("Press Enter to continue...")
            continue


if __name__ == "__main__":
    main()
