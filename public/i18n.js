// i18n.js

i18next.init({
    lng: 'fr', // langue par défaut
    resources: {
        en: {
            translation: {
                "title": "Exexe System",
                "backButton": "Back to Index",
                "saveButton": "Save",
                "codePageTitle": "Code Examples" ,
                "screenTitle": "Screen",
                "resultsDisplay": "The results will be displayed here",
                "goal_with_script": "Goal with predefined script",
                "execute": "Execute",
                "display_trace": "Display Trace",
                "display_examples": "Display Prolog Examples",
                "display_script": "Display Exexe Script",
                "trace_display": "The trace is displayed here:",
                "script_construction": "Script Construction",
                "initialization": "Initialization",
                "enter_initial_goal": "Enter an initial goal",
                "send": "Send",
                "goal_choice": "Goal Choice",
                "rule": "Rule",
                "structural_induction_variable": "Choice of structural induction variable",
                "confirmation": "Confirmation",
                "combinedScriptPageTitle": "Combined Script",
                "scriptTitle": "Script",
                "codeExamplesTitle": "Prolog Code Examples",
                "codePageTitle": "Code Examples Page",
               "enterFileNamePrompt": "Enter file name:"
            }
        },
        fr: {
            translation: {
                "title": "Le système Exexe",
                "backButton": "Revenir à l'index",
                "saveButton": "Enregistrer",
                "screenTitle": "Écran",
                "resultsDisplay": "Les résultats seront affichés ici",
                "goal_with_script": "But avec script prédéfini",
                "execute": "Exécuter",
                "display_trace": "Afficher la trace",
                "display_examples": "Afficher les exemples Prolog",
                "display_script": "Afficher le script Exexe",
                "trace_display": "La trace s'affiche ici :",
                "script_construction": "Construction de script",
                "initialization": "Initialisation",
                "enter_initial_goal": "Entrez un but initial",
                "send": "Envoyer",
                "goal_choice": "Choix du but",
                "rule": "Règle",
                "structural_induction_variable": "Choix de la variable d'induction structurelle",
                "confirmation": "Confirmation",
                "combinedScriptPageTitle": "Script Combiné",
                "scriptTitle": "Script",
                "codeExamplesTitle": "Exemples de Code Prolog",
                "codePageTitle": "Page des Exemples de Code",
                "enterFileNamePrompt": "Saisissez le nom du fichier :" 
            }
        }
    },
    fallbackLng: 'fr',
    debug: true
}, function(err, t) {
    // Initialisation terminée
    updateContent();
});

function updateContent() {
    document.querySelectorAll('[data-i18n]').forEach(function(element) {
        var key = element.getAttribute('data-i18n');
        element.innerHTML = i18next.t(key);
    });
}

function changeLanguage(lng) {
    i18next.changeLanguage(lng, function(err, t) {
        if (err) return console.log('something went wrong loading', err);
        updateContent();
    });
}
