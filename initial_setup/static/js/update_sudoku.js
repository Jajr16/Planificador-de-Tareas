$(document).ready(function(){
    $('#Res_Sudoku').on('input', '.sudoku-cell', function(){
        try {
            var sudoku_values = [];
            $('.sudoku-cell').each(function () {
                var valor = $(this).text().trim() || 0;
                sudoku_values.push(valor);
            });

            // Organiza los valores en una matriz 9x9
            var sudoku_matrix = [];
            for (var i = 0; i < 9; i++) {
                var row_values = sudoku_values.slice(i * 9, (i + 1) * 9);
                sudoku_matrix.push(row_values);
            }

            // Puedes imprimir la matriz en la consola para verificar
            console.log("Matriz 9x9:", sudoku_matrix);

            // Obtén el token CSRF del documento
            var csrftoken = getCookie('csrftoken');

            $.ajax({
                url: 'UpSudo',
                type: 'POST',
                contentType: 'application/json; charset=utf-8',
                data: JSON.stringify({ 'sudoku_values': sudoku_matrix }),
                headers: { 'X-CSRFToken': csrftoken },  // Agregar el token CSRF en la cabecera
                success: function (response) {
                    console.log("Respuesta AJAX:", response);
                    if (response.success === 0)
                        alert(response.message);
                },
                error: function (xhr, status, error) {
                    console.error("AJAX request failed:", status, error);
                }
            });
        } catch (error) {
            console.error("An error occurred:", error);
        }
    });
});

// Función para obtener el valor del token CSRF de la cookie
function getCookie(name) {
    var cookieValue = null;
    if (document.cookie && document.cookie !== '') {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = cookies[i].trim();
            // Buscar el nombre del cookie con el formato 'csrftoken='
            if (cookie.substring(0, name.length + 1) === (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}
