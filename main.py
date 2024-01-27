import sys
import time
sys.path.insert(0, './lib')

from rich.console import Console
from rich.markdown import Markdown
from rich.panel import Panel

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




def main():
    while True:
        choice = None
        while choice is None:
            choice = create_menu(options)

        if choice == 0:
         time.sleep(5)
        if choice == 1:
            time.sleep(5)

        if choice == 2:
          time.sleep(5)

        if choice == 3:
            
            time.sleep(5)    

if __name__ == "__main__":
    main()
    
